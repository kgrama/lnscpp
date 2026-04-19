// pbf.cpp — Parameterized Bounded Format (Table-Free Edition) in C++
// Port of pbf.py: pure-integer CF-based arithmetic, no lookup tables.
//
// Storage: a PBF code is an unsigned integer in [0, (1<<n_bits)-1].
//   code == 0            → zero
//   code == max_code     → infinity
//   code <  mid          → negative values (larger mag as code decreases)
//   code >= mid          → positive values (larger mag as code increases)
//
// Internally, CF evaluation uses Q32 fixed-point (Q = 2^32).

#pragma once

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>

// ── Fixed-point Q32 constants ────────────────────────────────────────────────

static const int64_t PBF_Q   = (int64_t)1 << 32;
static const int64_t PBF_LN2 = 2977044472LL;   // round(ln(2) * 2^32)

// ── Format descriptor ────────────────────────────────────────────────────────

/// Parameterized Bounded Format descriptor. Created once per format.
struct pbf_t {
    int n_bits;          // bit width (4, 8, 16, ...)
    double v_min;        // smallest representable magnitude
    double v_max;        // largest representable magnitude
    uint32_t max_code;   // (1 << n_bits) - 1, used as +inf sentinel
    uint32_t mid;        // 1 << (n_bits - 1), sign boundary
    int32_t  n_levels;   // mid - 1, distinct magnitude levels per sign
    int64_t  scale;      // Q-format: ln_range / n_levels (floor)
    int64_t  offset;     // Q-format: round(ln(v_min) / scale), used in mul/div
};

// ═════════════════════════════════════════════════════════════════════════════
// INTEGER HELPERS
// ═════════════════════════════════════════════════════════════════════════════

/// Count leading zeros of x in an n_bits-wide field. Returns n_bits if x==0.
static inline int pbf_clz(uint32_t x, int n_bits) {
    if (x == 0) return n_bits;
    int len = 0;
    while (x) { len++; x >>= 1; }
    return n_bits - len;
}

/// Python floor-division: floor(a/b). Handles negative operands correctly.
static inline int64_t pbf_floordiv(int64_t a, int64_t b) {
    int64_t q = a / b;
    int64_t r = a - q * b;
    if (r != 0 && ((r < 0) != (b < 0))) q -= 1;
    return q;
}

/// Kuttaka division: round a/b to nearest, ties toward +infinity.
/// Matches Python `(a + b//2) // b` for b>0 and `(a - b//2) // b` for b<0.
static inline int64_t pbf_kdiv(int64_t a, int64_t b) {
    if (b == 0) return a;
    int64_t b_half = pbf_floordiv(b, 2);
    int64_t num = (b > 0) ? (a + b_half) : (a - b_half);
    return pbf_floordiv(num, b);
}

/// kdiv(a*b, c), computed in __int128 to avoid overflow.
static inline int64_t pbf_muldiv(int64_t a, int64_t b, int64_t c) {
    if (c == 0) return a;
    __int128 prod  = (__int128)a * (__int128)b;
    __int128 chalf = (c > 0) ? (__int128)(c / 2) : -(__int128)((-c) / 2);
    __int128 num   = (c > 0) ? (prod + chalf) : (prod - chalf);
    __int128 q = num / (__int128)c;
    __int128 r = num - q * (__int128)c;
    if (r != 0 && ((r < 0) != (c < 0))) q -= 1;
    return (int64_t)q;
}

// ═════════════════════════════════════════════════════════════════════════════
// CONTINUED FRACTIONS
// ═════════════════════════════════════════════════════════════════════════════

/// ln(1 + x/Q) * Q via an 8-term continued fraction (Euler).
static inline int64_t pbf_ln1p_cf(int64_t x) {
    if (x == 0)        return 0;
    if (x <= -PBF_Q)   return -PBF_Q * 100;
    int64_t h = 8 * PBF_Q;
    h = 7 * PBF_Q + pbf_muldiv(16 * x, PBF_Q, h);
    h = 6 * PBF_Q + pbf_muldiv( 9 * x, PBF_Q, h);
    h = 5 * PBF_Q + pbf_muldiv( 9 * x, PBF_Q, h);
    h = 4 * PBF_Q + pbf_muldiv( 4 * x, PBF_Q, h);
    h = 3 * PBF_Q + pbf_muldiv( 4 * x, PBF_Q, h);
    h = 2 * PBF_Q + pbf_muldiv(     x, PBF_Q, h);
    h = 1 * PBF_Q + pbf_muldiv(     x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

/// ln(x/Q) * Q with range reduction. Accepts __int128 to handle wide inputs.
static inline int64_t pbf_ln_cf_big(__int128 x) {
    if (x <= 0)      return -PBF_Q * 100;
    if (x == PBF_Q)  return 0;
    int64_t k = 0;
    __int128 m = x;
    __int128 Q128 = (__int128)PBF_Q;
    while (m >= 2 * Q128) { m >>= 1; k++; }
    while (m < Q128 / 2)  { m <<= 1; k--; }
    int64_t m64 = (int64_t)m;
    return k * PBF_LN2 + pbf_ln1p_cf(m64 - PBF_Q);
}

/// ln(x/Q) * Q for int64 inputs (cannot represent very-large ratios).
static inline int64_t pbf_ln_cf(int64_t x) { return pbf_ln_cf_big((__int128)x); }

/// (exp(x/Q) - 1) * Q via a 7-term continued fraction.
static inline int64_t pbf_expm1_cf(int64_t x) {
    if (x == 0)           return 0;
    if (x >  20 * PBF_Q)  return PBF_Q * (1 << 20);
    if (x < -20 * PBF_Q)  return -PBF_Q;
    int64_t h = 2 * PBF_Q;
    h = 7 * PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2 * PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 5 * PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2 * PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 3 * PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2 * PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 1 * PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

/// exp(x/Q) * Q with ln2 range reduction.
static inline int64_t pbf_exp_cf(int64_t x) {
    if (x == 0) return PBF_Q;
    int64_t k  = pbf_kdiv(x, PBF_LN2);
    int64_t r  = x - k * PBF_LN2;
    int64_t er = PBF_Q + pbf_expm1_cf(r);
    if (k >= 0) {
        if (k >= 62) return PBF_Q * (1LL << 20);
        return er << k;
    }
    int64_t nk = -k;
    if (nk >= 63) return 0;
    return er >> nk;
}

/// ln(1 + exp(d/Q)) * Q — softplus in Q format.
static inline int64_t pbf_softplus_cf(int64_t d) {
    if (d >  10 * PBF_Q) return d;
    if (d < -10 * PBF_Q) return pbf_exp_cf(d);
    if (d > 0) {
        int64_t e = pbf_exp_cf(-d);
        return d + pbf_ln1p_cf(e);
    }
    int64_t e = pbf_exp_cf(d);
    return pbf_ln1p_cf(e);
}

/// ln(1 - exp(d/Q)) * Q — softminus, d < 0.
static inline int64_t pbf_softminus_cf(int64_t d) {
    if (d >= 0)          return -100 * PBF_Q;
    if (d < -10 * PBF_Q) return -pbf_exp_cf(d);
    int64_t e = pbf_exp_cf(d);
    return pbf_ln1p_cf(-e);
}

// ═════════════════════════════════════════════════════════════════════════════
// FORMAT INITIALIZATION
// ═════════════════════════════════════════════════════════════════════════════

/// Initialize a PBF descriptor with the given bit width and magnitude range.
static void pbf_init(pbf_t* p, int n_bits, double v_min, double v_max) {
    p->n_bits   = n_bits;
    p->v_min    = v_min;
    p->v_max    = v_max;
    p->max_code = ((uint32_t)1 << n_bits) - 1u;
    p->mid      = (uint32_t)1 << (n_bits - 1);
    p->n_levels = (int32_t)p->mid - 1;

    // Init constants are computed in double then converted to Q32. The CF
    // path is lossy here: representing tiny v_min in Q32 truncates to a few
    // significant digits (1e-8*Q rounds to 42, ≈2% off), and pbf_ln_cf_big
    // accumulates O(k) LSB error from k * PBF_LN2 at large arguments.
    // Runtime arithmetic still uses the CF — only the one-shot init constants
    // go through libm.
    double ln_vmin_d = std::log(v_min);
    double ln_vmax_d = std::log(v_max);
    int64_t ln_range = (int64_t)std::llround((ln_vmax_d - ln_vmin_d) * (double)PBF_Q);
    int64_t ln_vmin  = (int64_t)std::llround(ln_vmin_d * (double)PBF_Q);

    // scale = ln_range // n_levels (Python floor division)
    p->scale = pbf_floordiv(ln_range, (int64_t)p->n_levels);

    // offset = round(ln(v_min) / scale), used by mul/div
    p->offset = (p->scale != 0) ? pbf_kdiv(ln_vmin, p->scale) : 0;
}

/// Preset: PBF4 with v_min=0.1, v_max=10.
static inline void pbf4_init(pbf_t* p)  { pbf_init(p, 4,  0.1,  10.0);  }

/// Preset: PBF8 with v_min=1e-4, v_max=1e4.
static inline void pbf8_init(pbf_t* p)  { pbf_init(p, 8,  1e-4, 1e4);   }

/// Preset: PBF16 with v_min=1e-8, v_max=1e8.
static inline void pbf16_init(pbf_t* p) { pbf_init(p, 16, 1e-8, 1e8);   }

// ═════════════════════════════════════════════════════════════════════════════
// ENCODE / DECODE
// ═════════════════════════════════════════════════════════════════════════════

/// Convert a floating-point value to a PBF code.
static uint32_t pbf_encode(const pbf_t* p, double value) {
    if (value == 0.0)  return 0;
    if (std::isinf(value)) return p->max_code;

    int sign = (value > 0.0) ? 1 : -1;
    double mag = fabs(value);
    if (mag < p->v_min) mag = p->v_min;
    if (mag > p->v_max) mag = p->v_max;

    __int128 ratio_q = (__int128)((mag / p->v_min) * (double)PBF_Q);
    int64_t  ln_ratio = pbf_ln_cf_big(ratio_q);

    int64_t level = pbf_kdiv(ln_ratio, p->scale);
    if (level < 0) level = 0;
    if (level > (int64_t)p->n_levels - 1) level = p->n_levels - 1;

    if (sign > 0) return p->mid + (uint32_t)level;
    return p->mid - 1u - (uint32_t)level;
}

/// Convert a PBF code back to floating-point.
/// Uses double-precision exp at the float boundary to avoid int64 overflow
/// for wide formats (e.g. PBF16 where level*scale can exceed 2^40).
static double pbf_decode(const pbf_t* p, uint32_t code) {
    if (code == 0)             return 0.0;
    if (code == p->max_code)   return INFINITY;

    double sign;
    int64_t level;
    if (code >= p->mid) { level = (int64_t)code - (int64_t)p->mid; sign = +1.0; }
    else                { level = (int64_t)p->mid - 1 - (int64_t)code; sign = -1.0; }

    double ln_ratio = (double)(level * p->scale) / (double)PBF_Q;
    return sign * p->v_min * std::exp(ln_ratio);
}

// ═════════════════════════════════════════════════════════════════════════════
// ARITHMETIC
// ═════════════════════════════════════════════════════════════════════════════

/// Negation: ones'-complement within the format range.
static inline uint32_t pbf_neg(const pbf_t* p, uint32_t a) {
    if (a == 0 || a == p->max_code) return a;
    return p->max_code - a;
}

/// Extract magnitude level from a code.
static inline int64_t pbf_mag(const pbf_t* p, uint32_t a, int a_neg) {
    if (a_neg) return (int64_t)p->mid - 1 - (int64_t)a;
    return (int64_t)a - (int64_t)p->mid;
}

/// Compose a code from sign and magnitude level.
static inline uint32_t pbf_compose(const pbf_t* p, int neg, int64_t mag) {
    if (neg) return (uint32_t)((int64_t)p->mid - 1 - mag);
    return (uint32_t)((int64_t)p->mid + mag);
}

/// Multiplication: add magnitudes in log domain, apply v_min offset.
static uint32_t pbf_mul(const pbf_t* p, uint32_t a, uint32_t b) {
    if (a == 0 || b == 0) return 0;
    if (a == p->max_code || b == p->max_code) return p->max_code;

    int a_neg = (a < p->mid);
    int b_neg = (b < p->mid);
    int64_t a_mag = pbf_mag(p, a, a_neg);
    int64_t b_mag = pbf_mag(p, b, b_neg);

    int64_t c_mag = a_mag + b_mag + p->offset;
    if (c_mag < 0) return 0;
    if (c_mag >= (int64_t)p->n_levels) return p->max_code;

    return pbf_compose(p, a_neg != b_neg, c_mag);
}

/// Division: subtract magnitudes in log domain, remove v_min offset.
static uint32_t pbf_div(const pbf_t* p, uint32_t a, uint32_t b) {
    if (a == 0) return 0;
    if (b == 0 || a == p->max_code) return p->max_code;
    if (b == p->max_code) return 0;

    int a_neg = (a < p->mid);
    int b_neg = (b < p->mid);
    int64_t a_mag = pbf_mag(p, a, a_neg);
    int64_t b_mag = pbf_mag(p, b, b_neg);

    int64_t c_mag = a_mag - b_mag - p->offset;
    if (c_mag < 0) return 0;
    if (c_mag >= (int64_t)p->n_levels) return p->max_code;

    return pbf_compose(p, a_neg != b_neg, c_mag);
}

/// Same-sign addition: magnitude combine via softplus CF.
static uint32_t pbf_add_same_sign(const pbf_t* p, uint32_t a, uint32_t b, int is_neg) {
    int64_t a_mag = pbf_mag(p, a, is_neg);
    int64_t b_mag = pbf_mag(p, b, is_neg);
    if (a_mag < b_mag) { int64_t t = a_mag; a_mag = b_mag; b_mag = t; }

    int64_t d = (b_mag - a_mag) * p->scale;
    int64_t correction = pbf_softplus_cf(d);
    int64_t c_levels   = pbf_kdiv(correction, p->scale);

    int64_t c_mag = a_mag + c_levels;
    if (c_mag < 0) c_mag = 0;
    if (c_mag > (int64_t)p->n_levels - 1) c_mag = p->n_levels - 1;
    return pbf_compose(p, is_neg, c_mag);
}

/// Opposite-sign addition: magnitude subtract via softminus CF.
static uint32_t pbf_add_diff_sign(const pbf_t* p, uint32_t a, uint32_t b,
                                  int a_neg, int b_neg) {
    int64_t a_mag = pbf_mag(p, a, a_neg);
    int64_t b_mag = pbf_mag(p, b, b_neg);
    if (a_mag == b_mag) return 0;

    int64_t big_mag, small_mag;
    int result_neg;
    if (a_mag > b_mag) { big_mag = a_mag; small_mag = b_mag; result_neg = a_neg; }
    else               { big_mag = b_mag; small_mag = a_mag; result_neg = b_neg; }

    int64_t d = (small_mag - big_mag) * p->scale;     // d < 0
    int64_t correction = pbf_softminus_cf(d);
    int64_t c_levels   = pbf_kdiv(correction, p->scale);

    int64_t c_mag = big_mag + c_levels;
    if (c_mag < 0) c_mag = 0;
    if (c_mag > (int64_t)p->n_levels - 1) c_mag = p->n_levels - 1;
    return pbf_compose(p, result_neg, c_mag);
}

/// Addition: dispatches on sign match.
static uint32_t pbf_add(const pbf_t* p, uint32_t a, uint32_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == p->max_code || b == p->max_code) return p->max_code;

    int a_neg = (a < p->mid);
    int b_neg = (b < p->mid);
    if (a_neg == b_neg) return pbf_add_same_sign(p, a, b, a_neg);
    return pbf_add_diff_sign(p, a, b, a_neg, b_neg);
}

/// Subtraction via add-negate.
static inline uint32_t pbf_sub(const pbf_t* p, uint32_t a, uint32_t b) {
    return pbf_add(p, a, pbf_neg(p, b));
}

/// Similarity: CLZ of XOR within the format width.
static inline int pbf_similarity(const pbf_t* p, uint32_t a, uint32_t b) {
    return pbf_clz(a ^ b, p->n_bits);
}
