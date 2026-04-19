// pbf_taper.cpp — Tapered PBF: dual-density log encoding
//
// Concentrates codes near unity (log magnitude 0, value 1.0) by partitioning
// the magnitude code space into 4 tapers with non-uniform log step sizes.
// Inner tapers (around the unity point) have a fine step; outer tapers
// (near v_min and v_max) have a coarser step. The result: better SNR for
// operands near 1.0, gracefully degrading for extreme magnitudes — matching
// the typical operand distribution in NN workloads, signal processing, etc.
//
// Layout (W16, max_mag = 32766, step ratios T0:T1:T2:T3 = 8:1:1:8):
//
//   T0 (coarse, sub-unity tail):  codes [0,     8191],  step = 8s
//   T1 (fine,   sub-unity):       codes [8192,  16382], step = s
//   T2 (fine,   super-unity):     codes [16383, 24573], step = s
//   T3 (coarse, super-unity tail):codes [24574, 32766], step = 8s
//
// where s is chosen so total log range = 2 * ln(v_max).
//
// The boundaries (8191, 16382, 24573, 32766) match P32_Taper.v exactly, so
// the structural soundness proofs (totality, disjointness, monotonicity,
// reconstructability) carry over unchanged. Only the per-taper step sizes
// — an interpretation-layer choice — differ from the proof's defaults.
//
// The format is symmetric: value 1.0 sits at the boundary between T1 and T2,
// the smallest representable positive magnitude is 1/v_max, the largest
// is v_max. Sign handled via the standard SBP mirror (negation = ones'
// complement of the code).
//
// Arithmetic: codes are decoded to a Q32 log magnitude, the operation runs
// in Q32 (integer add for mul/div, sb/db CF for add/sub), then re-encoded.
// This adds two encode/decode passes per op vs flat PBF, but the encoding
// gain at unity dominates the round-trip cost.

#include "pbf.cpp"

#define PBF_TAPER_N 4

typedef struct {
    int      n_bits;
    uint32_t mid;
    uint32_t max_code;
    int32_t  max_mag;                    // = mid - 2 (largest valid magnitude code)

    int32_t  taper_end[PBF_TAPER_N];     // inclusive end of each taper (in mag code units)
    int64_t  taper_step[PBF_TAPER_N];    // Q32 log step within each taper
    int64_t  taper_log[PBF_TAPER_N + 1]; // Q32 cumulative log magnitude at the
                                         // boundary preceding taper k (taper_log[0]
                                         // = log at mag code 0 = -ln(v_max))

    double   v_max;
} pbf_taper_t;

// ─────────────────────────────────────────────────────────────────────────────
// FORMAT INITIALIZATION
// ─────────────────────────────────────────────────────────────────────────────

/// Initialize a tapered PBF descriptor. Currently W16-only (matches the
/// boundaries in proofs/P32_Taper.v); other widths could be supported by
/// scaling per_inner = max_mag / 4.
static void pbf_taper_init(pbf_taper_t* p, int n_bits, double v_max) {
    p->n_bits   = n_bits;
    p->mid      = (uint32_t)1 << (n_bits - 1);
    p->max_code = ((uint32_t)1 << n_bits) - 1u;
    p->max_mag  = (int32_t)p->mid - 2;
    p->v_max    = v_max;

    if (n_bits == 16) {
        // Match P32_Taper.v constants exactly so the proof applies verbatim.
        p->taper_end[0] = 8191;
        p->taper_end[1] = 16382;
        p->taper_end[2] = 24573;
        p->taper_end[3] = 32766;
    } else {
        // Generic equal-quarters partition with the residue absorbed by T3.
        int32_t per = p->max_mag / 4;
        p->taper_end[0] = per - 1;
        p->taper_end[1] = 2 * per - 1;
        p->taper_end[2] = 3 * per - 1;
        p->taper_end[3] = p->max_mag;
    }

    // Step ratios: outer tapers 8x coarser than inner. This places ~88% of
    // the codes in the outer tapers (8x cheaper to span) while the inner
    // tapers cover the high-precision region near unity.
    static const int32_t RATIOS[PBF_TAPER_N] = { 8, 1, 1, 8 };

    // Solve for the base step s such that the total log span equals 2*ln(v_max).
    int32_t prev_end = -1;
    int64_t weighted = 0;
    for (int k = 0; k < PBF_TAPER_N; k++) {
        int32_t codes = p->taper_end[k] - prev_end;
        weighted += (int64_t)codes * RATIOS[k];
        prev_end = p->taper_end[k];
    }
    double base_step = 2.0 * std::log(v_max) / (double)weighted;

    // Build per-taper Q32 step and cumulative log boundaries. The format is
    // symmetric in log space: taper_log[0] = -ln(v_max), taper_log[4] = +ln(v_max).
    p->taper_log[0] = -(int64_t)std::llround(std::log(v_max) * (double)PBF_Q);
    prev_end = -1;
    for (int k = 0; k < PBF_TAPER_N; k++) {
        p->taper_step[k] = (int64_t)std::llround(base_step * RATIOS[k] * (double)PBF_Q);
        int32_t codes = p->taper_end[k] - prev_end;
        p->taper_log[k + 1] = p->taper_log[k] + p->taper_step[k] * (int64_t)codes;
        prev_end = p->taper_end[k];
    }
}

/// Preset: tapered PBF8 with v_max = 1e4 (range [1e-4, 1e4]).
static inline void pbf_taper8_init (pbf_taper_t* p) { pbf_taper_init(p,  8, 1e4); }

/// Preset: tapered PBF16 with v_max = 1e8 (range [1e-8, 1e8]).
static inline void pbf_taper16_init(pbf_taper_t* p) { pbf_taper_init(p, 16, 1e8); }

// ─────────────────────────────────────────────────────────────────────────────
// MAGNITUDE CODE ↔ Q32 LOG MAGNITUDE
// ─────────────────────────────────────────────────────────────────────────────

/// Find the taper index containing a magnitude code. Linear scan over 4
/// boundaries — branches predict well and beat binary search at this size.
static inline int pbf_taper_find_by_mag(const pbf_taper_t* p, int32_t mag) {
    if (mag <= p->taper_end[0]) return 0;
    if (mag <= p->taper_end[1]) return 1;
    if (mag <= p->taper_end[2]) return 2;
    return 3;
}

/// Find the taper index containing a Q32 log magnitude.
static inline int pbf_taper_find_by_log(const pbf_taper_t* p, int64_t log_q32) {
    if (log_q32 <= p->taper_log[1]) return 0;
    if (log_q32 <= p->taper_log[2]) return 1;
    if (log_q32 <= p->taper_log[3]) return 2;
    return 3;
}

/// Decode a magnitude code to its signed Q32 log magnitude.
static inline int64_t pbf_taper_mag_to_log(const pbf_taper_t* p, int32_t mag) {
    int     k    = pbf_taper_find_by_mag(p, mag);
    int32_t base = (k == 0) ? 0 : (p->taper_end[k - 1] + 1);
    int32_t pos  = mag - base;
    return p->taper_log[k] + (int64_t)pos * p->taper_step[k];
}

/// Encode a signed Q32 log magnitude to the nearest magnitude code.
/// Clamps to [0, max_mag] on overflow.
static inline int32_t pbf_taper_log_to_mag(const pbf_taper_t* p, int64_t log_q32) {
    if (log_q32 <= p->taper_log[0]) return 0;
    if (log_q32 >= p->taper_log[4]) return p->max_mag;

    int     k    = pbf_taper_find_by_log(p, log_q32);
    int32_t base = (k == 0) ? 0 : (p->taper_end[k - 1] + 1);
    int64_t pos  = pbf_kdiv(log_q32 - p->taper_log[k], p->taper_step[k]);
    int32_t mag  = base + (int32_t)pos;
    if (mag < 0)             return 0;
    if (mag > p->max_mag)    return p->max_mag;
    return mag;
}

// ─────────────────────────────────────────────────────────────────────────────
// ENCODE / DECODE
// ─────────────────────────────────────────────────────────────────────────────

static uint32_t pbf_taper_encode(const pbf_taper_t* p, double value) {
    if (value == 0.0)      return 0;
    if (std::isinf(value)) return p->max_code;

    int    sign  = (value > 0.0) ? 1 : -1;
    double abs_v = std::fabs(value);
    if (abs_v > p->v_max)        abs_v = p->v_max;
    if (abs_v < 1.0 / p->v_max)  abs_v = 1.0 / p->v_max;

    int64_t log_q32 = (int64_t)std::llround(std::log(abs_v) * (double)PBF_Q);
    int32_t mag     = pbf_taper_log_to_mag(p, log_q32);

    if (sign > 0) return p->mid + (uint32_t)mag;
    return p->mid - 1u - (uint32_t)mag;
}

static double pbf_taper_decode(const pbf_taper_t* p, uint32_t code) {
    if (code == 0)            return 0.0;
    if (code == p->max_code)  return INFINITY;

    double  sign;
    int32_t mag;
    if (code >= p->mid) { mag = (int32_t)code - (int32_t)p->mid;       sign = +1.0; }
    else                { mag = (int32_t)p->mid - 1 - (int32_t)code;   sign = -1.0; }

    int64_t log_q32 = pbf_taper_mag_to_log(p, mag);
    return sign * std::exp((double)log_q32 / (double)PBF_Q);
}

// ─────────────────────────────────────────────────────────────────────────────
// ARITHMETIC
// ─────────────────────────────────────────────────────────────────────────────

static inline uint32_t pbf_taper_neg(const pbf_taper_t* p, uint32_t a) {
    if (a == 0 || a == p->max_code) return a;
    return p->max_code - a;
}

/// Extract (sign, signed Q32 log magnitude) from a code. Returns true if
/// the code is finite and non-zero (i.e., the log magnitude is meaningful).
static inline bool pbf_taper_unpack(const pbf_taper_t* p, uint32_t code,
                                    int* sign_neg, int64_t* log_q32) {
    if (code == 0 || code == p->max_code) return false;
    *sign_neg = (code < p->mid);
    int32_t mag = (*sign_neg) ? ((int32_t)p->mid - 1 - (int32_t)code)
                              : ((int32_t)code - (int32_t)p->mid);
    *log_q32 = pbf_taper_mag_to_log(p, mag);
    return true;
}

/// Compose a code from a sign flag and a signed Q32 log magnitude.
static inline uint32_t pbf_taper_pack(const pbf_taper_t* p, int sign_neg, int64_t log_q32) {
    int32_t mag = pbf_taper_log_to_mag(p, log_q32);
    if (sign_neg) return p->mid - 1u - (uint32_t)mag;
    return p->mid + (uint32_t)mag;
}

static uint32_t pbf_taper_mul(const pbf_taper_t* p, uint32_t a, uint32_t b) {
    if (a == 0 || b == 0) return 0;
    if (a == p->max_code || b == p->max_code) return p->max_code;

    int     a_neg, b_neg;
    int64_t la, lb;
    pbf_taper_unpack(p, a, &a_neg, &la);
    pbf_taper_unpack(p, b, &b_neg, &lb);
    return pbf_taper_pack(p, a_neg ^ b_neg, la + lb);
}

static uint32_t pbf_taper_div(const pbf_taper_t* p, uint32_t a, uint32_t b) {
    if (a == 0)              return 0;
    if (b == 0 || a == p->max_code) return p->max_code;
    if (b == p->max_code)    return 0;

    int     a_neg, b_neg;
    int64_t la, lb;
    pbf_taper_unpack(p, a, &a_neg, &la);
    pbf_taper_unpack(p, b, &b_neg, &lb);
    return pbf_taper_pack(p, a_neg ^ b_neg, la - lb);
}

static uint32_t pbf_taper_add(const pbf_taper_t* p, uint32_t a, uint32_t b) {
    if (a == 0)            return b;
    if (b == 0)            return a;
    if (a == p->max_code || b == p->max_code) return p->max_code;

    int     a_neg, b_neg;
    int64_t la, lb;
    pbf_taper_unpack(p, a, &a_neg, &la);
    pbf_taper_unpack(p, b, &b_neg, &lb);

    // Order so |a| >= |b| (la >= lb). Sign of result follows the larger operand
    // for opposite-sign add; for same-sign it follows either.
    int big_neg, small_neg;
    int64_t big_log, small_log;
    if (la >= lb) { big_log = la; small_log = lb; big_neg = a_neg; small_neg = b_neg; }
    else          { big_log = lb; small_log = la; big_neg = b_neg; small_neg = a_neg; }

    int64_t d = small_log - big_log;   // <= 0

    int64_t log_c;
    int     c_neg;
    if (big_neg == small_neg) {
        // |a|+|b|: log_c = log|a| + log(1 + exp(log|b|-log|a|))
        log_c = big_log + pbf_softplus_cf(d);
        c_neg = big_neg;
    } else {
        // |a|-|b|: log_c = log|a| + log(1 - exp(log|b|-log|a|))
        if (d == 0) return 0;          // exact cancellation
        log_c = big_log + pbf_softminus_cf(d);
        c_neg = big_neg;
    }
    return pbf_taper_pack(p, c_neg, log_c);
}

static inline uint32_t pbf_taper_sub(const pbf_taper_t* p, uint32_t a, uint32_t b) {
    return pbf_taper_add(p, a, pbf_taper_neg(p, b));
}
