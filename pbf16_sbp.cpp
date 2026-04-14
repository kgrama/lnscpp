// pbf16_sbp.cpp — Full SBP (Sign-By-Position) LNS implementation
//
// Patent pending GB 2602876.1, snums (2026)
//
// LICENSE: Research and evaluation use only.
// Commercial use is explicitly prohibited without a separate written
// license from the copyright holder. See LICENSE for full terms.
//
// FULL SBP ENCODING (snums Claim 1):
//   - NOT(x) = -x (single gate negation)
//   - Sign determined by position relative to mirror axis
//   - CLZ-adaptive CF depth (snums Claim 12)
//
// API compatible with xlns16.cpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#ifdef _WIN32
typedef unsigned __int16 xlns16;
typedef __int16          xlns16_signed;
#else
#include <sys/types.h>
typedef u_int16_t xlns16;
typedef int16_t   xlns16_signed;
#endif

// ══════════════════════════════════════════════════════════════════════════════
// SBP ENCODING (snums Claim 1)
// ══════════════════════════════════════════════════════════════════════════════
//
// 16-bit SBP layout with mirror axis at 0x7FFF.5:
//
//   Code 0x0000 = Zero (positive boundary)
//   Code 0x0001..0x7FFF = Negative values
//        0x0001 = smallest negative magnitude (closest to zero)
//        0x7FFF = largest negative magnitude (-Inf)
//   Code 0x8000..0xFFFE = Positive values  
//        0x8000 = largest positive magnitude (+Inf)
//        0xFFFE = smallest positive magnitude (closest to zero)
//   Code 0xFFFF = NaN (= ~0x0000)
//
// Mirror symmetry: ~x maps x to -x
//   ~0x0000 = 0xFFFF (zero ↔ NaN)
//   ~0x7FFF = 0x8000 (-Inf ↔ +Inf)  
//   ~0x0001 = 0xFFFE (smallest neg ↔ smallest pos)
//
// Log magnitude encoding:
//   For positive (code >= 0x8000): log2(|v|) * 128 = 0xBFFF - code + 1
//   For negative (code < 0x8000):  log2(|v|) * 128 = code - 1
//   Code 0x8000 + 0x3FFF = 0xBFFF encodes log2(1) = 0 → value = +1.0
//   Code 0x7FFF - 0x3FFF = 0x4000 encodes log2(1) = 0 → value = -1.0

#define SBP_ZERO        0x0000u
#define SBP_NAN         0xFFFFu
#define SBP_NEG_INF     0x7FFFu
#define SBP_POS_INF     0x8000u
#define SBP_AXIS        0x8000u   // Values >= this are positive
#define SBP_ONE         0xBFFFu   // +1.0: distance 0x3FFF from axis
#define SBP_NEG_ONE     0x4000u   // -1.0: ~SBP_ONE = 0x4000
#define SBP_SCALE       128       // log2 fractional bits

// Offset for log encoding: distance 0x3FFF = log2(1) = 0
#define SBP_LOG_OFFSET  0x3FFFu

// ── SBP core operations (snums Claim 1) ──────────────────────────────────────

// Negation: NOT (single gate) — THE key SBP property
#define sbp_neg(x)          ((xlns16)(~(x)))

// Sign extraction
#define sbp_is_zero(x)      ((x) == SBP_ZERO)
#define sbp_is_nan(x)       ((x) == SBP_NAN)
#define sbp_is_positive(x)  ((x) >= SBP_AXIS && (x) != SBP_NAN)
#define sbp_is_negative(x)  ((x) > SBP_ZERO && (x) < SBP_AXIS)

// Magnitude distance from mirror axis (unsigned)
// Positive: code - 0x8000, Negative: 0x7FFF - code
static inline uint16_t sbp_distance(xlns16 x) {
    if (x >= SBP_AXIS) return x - SBP_AXIS;
    return SBP_AXIS - 1 - x;
}

// Absolute value: map negative to positive via NOT
static inline xlns16 sbp_abs(xlns16 x) {
    return sbp_is_negative(x) ? sbp_neg(x) : x;
}

// Log magnitude (signed): distance - offset
// Returns log2(|v|) * 128 as signed value
static inline int16_t sbp_log_mag(xlns16 x) {
    return (int16_t)(sbp_distance(x)) - (int16_t)SBP_LOG_OFFSET;
}

// Construct SBP code from sign and log magnitude
static inline xlns16 sbp_from_sign_log(int negative, int16_t log_mag) {
    uint16_t dist = (uint16_t)((int16_t)SBP_LOG_OFFSET + log_mag);
    if (dist > 0x7FFF) dist = 0x7FFF;  // clamp to max
    if (negative) {
        return (xlns16)(SBP_AXIS - 1 - dist);
    } else {
        return (xlns16)(SBP_AXIS + dist);
    }
}

// ── CLZ for adaptive precision (snums Claim 12) ──────────────────────────────

static inline int sbp_clz16(uint16_t x) {
    if (x == 0) return 16;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clz((unsigned int)x) - 16;
#else
    int n = 0;
    if (x <= 0x00FF) { n += 8; x <<= 8; }
    if (x <= 0x0FFF) { n += 4; x <<= 4; }
    if (x <= 0x3FFF) { n += 2; x <<= 2; }
    if (x <= 0x7FFF) { n += 1; }
    return n;
#endif
}

// Magnitude delta for cancellation detection (snums Claim 8)
static inline uint16_t sbp_mag_delta(xlns16 a, xlns16 b) {
    uint16_t da = sbp_distance(a);
    uint16_t db = sbp_distance(b);
    return (da > db) ? (da - db) : (db - da);
}

// ══════════════════════════════════════════════════════════════════════════════
// LEGACY xlns16 API MAPPING
// ══════════════════════════════════════════════════════════════════════════════

#define xlns16_zero         SBP_ZERO
#define xlns16_one          SBP_ONE
#define xlns16_neg_one      SBP_NEG_ONE
#define xlns16_scale        SBP_SCALE
#define xlns16_logsignmask  SBP_LOG_OFFSET

// Legacy macros mapped to SBP
#define xlns16_neg(x)       sbp_neg(x)
#define xlns16_abs(x)       sbp_abs(x)
#define xlns16_sign(x)      (sbp_is_negative(x) ? 1 : 0)

// ══════════════════════════════════════════════════════════════════════════════
// PBF CONTINUED FRACTION ENGINE
// ══════════════════════════════════════════════════════════════════════════════

static const int64_t PBF_Q   = (int64_t)1 << 32;
static const int64_t PBF_LN2 = 2977044472LL;   // ln(2) * 2^32

// Kuttaka division: round-to-nearest
static inline int64_t pbf_kdiv(int64_t a, int64_t b) {
    if (b == 0) return a;
    if (b > 0) return (a + (a >= 0 ? b/2 : -(b/2))) / b;
    else       return (a + (a >= 0 ? -((-b)/2) : (-b)/2)) / b;
}

// Widening multiply-divide
static inline int64_t pbf_muldiv(int64_t a, int64_t b, int64_t c) {
    if (c == 0) return a;
    __int128 q = (__int128)a * b / c;
    __int128 r = (__int128)a * b - q * c;
    __int128 ac = c > 0 ? c : -c;
    __int128 ar = r > 0 ? r : -r;
    if (2*ar >= ac) q += (r*c > 0) ? 1 : -1;
    return (int64_t)q;
}

// ── CLZ-adaptive CF depth (snums Claim 12) ───────────────────────────────────

typedef enum {
    CF_DEPTH_8  = 8,
    CF_DEPTH_12 = 12,
    CF_DEPTH_16 = 16
} cf_depth_t;

static inline cf_depth_t sbp_select_cf_depth(uint16_t mag_delta) {
    int clz = sbp_clz16(mag_delta);
    if (clz < 8)  return CF_DEPTH_8;
    if (clz < 12) return CF_DEPTH_12;
    return CF_DEPTH_16;
}

// ── ln(1+x) with adaptive depth ──────────────────────────────────────────────

static inline int64_t pbf_ln1p_8(int64_t x) {
    if (x == 0)      return 0;
    if (x <= -PBF_Q) return -PBF_Q * 100;
    int64_t h = 8*PBF_Q;
    h = 7*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h = 6*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h = 5*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h = 4*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h = 3*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h = 2*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    h = 1*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_ln1p_12(int64_t x) {
    if (x == 0)      return 0;
    if (x <= -PBF_Q) return -PBF_Q * 100;
    int64_t h = 12*PBF_Q;
    h = 11*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);
    h = 10*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  9*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  8*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  7*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  6*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  5*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  4*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  3*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    h =  1*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_ln1p_16(int64_t x) {
    if (x == 0)      return 0;
    if (x <= -PBF_Q) return -PBF_Q * 100;
    int64_t h = 16*PBF_Q;
    h = 15*PBF_Q + pbf_muldiv(64*x, PBF_Q, h);
    h = 14*PBF_Q + pbf_muldiv(49*x, PBF_Q, h);
    h = 13*PBF_Q + pbf_muldiv(49*x, PBF_Q, h);
    h = 12*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);
    h = 11*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);
    h = 10*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  9*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  8*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  7*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  6*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  5*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  4*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  3*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    h =  1*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_ln1p(int64_t x) { return pbf_ln1p_8(x); }

static inline int64_t pbf_ln1p_adaptive(int64_t x, cf_depth_t depth) {
    switch (depth) {
        case CF_DEPTH_16: return pbf_ln1p_16(x);
        case CF_DEPTH_12: return pbf_ln1p_12(x);
        default:          return pbf_ln1p_8(x);
    }
}

// ── exp(x)-1 with adaptive depth ─────────────────────────────────────────────

static inline int64_t pbf_expm1_8(int64_t x) {
    if (x == 0)        return 0;
    if (x >  20*PBF_Q) return PBF_Q * (1 << 20);
    if (x < -20*PBF_Q) return -PBF_Q;
    int64_t h = 2*PBF_Q;
    h = 7*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 5*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 3*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h = 2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 1*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_expm1_12(int64_t x) {
    if (x == 0)        return 0;
    if (x >  20*PBF_Q) return PBF_Q * (1 << 20);
    if (x < -20*PBF_Q) return -PBF_Q;
    int64_t h = 2*PBF_Q;
    h = 11*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  9*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  7*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  5*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  3*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  1*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_expm1_16(int64_t x) {
    if (x == 0)        return 0;
    if (x >  20*PBF_Q) return PBF_Q * (1 << 20);
    if (x < -20*PBF_Q) return -PBF_Q;
    int64_t h = 2*PBF_Q;
    h = 15*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 13*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h = 11*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  9*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  7*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  5*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  3*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv( x, PBF_Q, h);
    h =  1*PBF_Q + pbf_muldiv(-x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

static inline int64_t pbf_expm1(int64_t x) { return pbf_expm1_8(x); }

static inline int64_t pbf_expm1_adaptive(int64_t x, cf_depth_t depth) {
    switch (depth) {
        case CF_DEPTH_16: return pbf_expm1_16(x);
        case CF_DEPTH_12: return pbf_expm1_12(x);
        default:          return pbf_expm1_8(x);
    }
}

// ── exp and ln with range reduction ──────────────────────────────────────────

static inline int64_t pbf_ln(int64_t x) {
    if (x <= 0)     return -PBF_Q * 100;
    if (x == PBF_Q) return 0;
    int64_t k = 0, m = x;
    while (m >= 2*PBF_Q) { m >>= 1; k++; }
    while (m <   PBF_Q/2) { m <<= 1; k--; }
    return k * PBF_LN2 + pbf_ln1p(m - PBF_Q);
}

static inline int64_t pbf_exp(int64_t x) {
    if (x == 0) return PBF_Q;
    int64_t k  = pbf_kdiv(x, PBF_LN2);
    int64_t r  = x - k * PBF_LN2;
    int64_t er = PBF_Q + pbf_expm1(r);
    return (k >= 0) ? (er << k) : (er >> (-k));
}

// ── softplus / softminus ─────────────────────────────────────────────────────

static inline int64_t pbf_softplus(int64_t d) {
    if (d >  10*PBF_Q) return d;
    if (d < -10*PBF_Q) return pbf_exp(d);
    return (d > 0) ? d + pbf_ln1p(pbf_exp(-d))
                   : pbf_ln1p(pbf_exp(d));
}

static inline int64_t pbf_softplus_adaptive(int64_t d, cf_depth_t depth) {
    if (d >  10*PBF_Q) return d;
    if (d < -10*PBF_Q) return pbf_exp(d);
    if (d > 0) {
        int64_t e = pbf_exp(-d);
        return d + pbf_ln1p_adaptive(e, depth);
    } else {
        int64_t e = pbf_exp(d);
        return pbf_ln1p_adaptive(e, depth);
    }
}

static inline int64_t pbf_softminus(int64_t d) {
    if (d >= 0)        return -100 * PBF_Q;
    if (d < -10*PBF_Q) return -pbf_exp(d);
    return pbf_ln1p(-pbf_exp(d));
}

static inline int64_t pbf_softminus_adaptive(int64_t d, cf_depth_t depth) {
    if (d >= 0)        return -100 * PBF_Q;
    if (d < -10*PBF_Q) return -pbf_exp(d);
    return pbf_ln1p_adaptive(-pbf_exp(d), depth);
}

// ══════════════════════════════════════════════════════════════════════════════
// SBP sb AND db (Gaussian logarithm replacement)
// ══════════════════════════════════════════════════════════════════════════════
//
// sb(z) = log2(1 + 2^(-z)) for same-sign addition
// db(z) = log2(1 - 2^(-z)) for opposite-sign (subtraction), z > 0

static inline int16_t sbp_sb(int16_t z) {
    // z is log magnitude difference (always >= 0 for sb)
    int64_t d  = pbf_kdiv((int64_t)z * PBF_LN2, SBP_SCALE);
    int64_t sp = pbf_softplus(-d);  // ln(1 + e^(-d))
    return (int16_t)pbf_kdiv(sp * SBP_SCALE, PBF_LN2);
}

static inline int16_t sbp_sb_adaptive(int16_t z, cf_depth_t depth) {
    int64_t d  = pbf_kdiv((int64_t)z * PBF_LN2, SBP_SCALE);
    int64_t sp = pbf_softplus_adaptive(-d, depth);
    return (int16_t)pbf_kdiv(sp * SBP_SCALE, PBF_LN2);
}

static inline int16_t sbp_db(int16_t z) {
    // z is log magnitude difference, must be > 0
    // db(z) = log2(1 - 2^(-z)) in scaled units
    if (z <= 0) return -32768;  // undefined
    int64_t d  = pbf_kdiv(-(int64_t)z * PBF_LN2, SBP_SCALE);
    int64_t sm = pbf_softminus(d);  // ln(1 - e^d) = ln(1 - 2^(-z))
    // Convert from ln domain to log2 domain: multiply by log2(e) = 1/ln(2)
    return (int16_t)pbf_kdiv(sm * SBP_SCALE, PBF_LN2);
}

static inline int16_t sbp_db_adaptive(int16_t z, cf_depth_t depth) {
    if (z <= 0) return -32768;
    int64_t d  = pbf_kdiv(-(int64_t)z * PBF_LN2, SBP_SCALE);
    int64_t sm = pbf_softminus_adaptive(d, depth);
    return (int16_t)pbf_kdiv(sm * SBP_SCALE, PBF_LN2);
}

// ══════════════════════════════════════════════════════════════════════════════
// ARITHMETIC OPERATIONS
// ══════════════════════════════════════════════════════════════════════════════

// Multiplication: add log magnitudes, XOR signs (via SBP: result sign from positions)
xlns16 xlns16_mul(xlns16 a, xlns16 b) {
    if (sbp_is_zero(a) || sbp_is_zero(b)) return SBP_ZERO;
    if (sbp_is_nan(a) || sbp_is_nan(b))   return SBP_NAN;
    
    int16_t log_a = sbp_log_mag(a);
    int16_t log_b = sbp_log_mag(b);
    int32_t log_sum = (int32_t)log_a + (int32_t)log_b;
    
    // Clamp to representable range
    if (log_sum > 0x3FFF)  log_sum = 0x3FFF;   // overflow to inf
    if (log_sum < -0x3FFF) log_sum = -0x3FFF;  // underflow to zero
    
    // Sign: negative if exactly one operand is negative
    int result_neg = sbp_is_negative(a) ^ sbp_is_negative(b);
    return sbp_from_sign_log(result_neg, (int16_t)log_sum);
}

// Division: subtract log magnitudes
xlns16 xlns16_div(xlns16 a, xlns16 b) {
    if (sbp_is_zero(b))                   return SBP_NAN;  // div by zero
    if (sbp_is_zero(a))                   return SBP_ZERO;
    if (sbp_is_nan(a) || sbp_is_nan(b))   return SBP_NAN;
    
    int16_t log_a = sbp_log_mag(a);
    int16_t log_b = sbp_log_mag(b);
    int32_t log_diff = (int32_t)log_a - (int32_t)log_b;
    
    if (log_diff > 0x3FFF)  log_diff = 0x3FFF;
    if (log_diff < -0x3FFF) log_diff = -0x3FFF;
    
    int result_neg = sbp_is_negative(a) ^ sbp_is_negative(b);
    return sbp_from_sign_log(result_neg, (int16_t)log_diff);
}

// Addition with CLZ-adaptive CF depth (snums Claim 12)
xlns16 xlns16_add(xlns16 a, xlns16 b) {
    // Handle special cases
    if (sbp_is_zero(a)) return b;
    if (sbp_is_zero(b)) return a;
    if (sbp_is_nan(a) || sbp_is_nan(b)) return SBP_NAN;
    
    // Get log magnitudes
    int16_t log_a = sbp_log_mag(a);
    int16_t log_b = sbp_log_mag(b);
    
    // Ensure |a| >= |b| (log_a >= log_b)
    if (log_a < log_b) {
        xlns16 t = a; a = b; b = t;
        int16_t tl = log_a; log_a = log_b; log_b = tl;
    }
    
    int16_t z = log_a - log_b;  // always >= 0
    
    // CLZ-adaptive depth selection
    cf_depth_t depth = sbp_select_cf_depth((uint16_t)z);
    
    int a_neg = sbp_is_negative(a);
    int b_neg = sbp_is_negative(b);
    
    if (a_neg == b_neg) {
        // Same sign: use sb (addition in magnitude)
        int16_t sb_val = sbp_sb_adaptive(z, depth);
        int16_t result_log = log_a + sb_val;
        if (result_log > 0x3FFF) result_log = 0x3FFF;
        return sbp_from_sign_log(a_neg, result_log);
    } else {
        // Opposite sign: use db (subtraction in magnitude)
        if (z == 0) return SBP_ZERO;  // exact cancellation
        int16_t db_val = sbp_db_adaptive(z, depth);
        int16_t result_log = log_a + db_val;
        // Result has sign of larger magnitude operand (a)
        return sbp_from_sign_log(a_neg, result_log);
    }
}

// Subtraction: add negated
#define xlns16_sub(a, b) xlns16_add((a), sbp_neg(b))

// ══════════════════════════════════════════════════════════════════════════════
// FLOAT ↔ SBP CONVERSION
// ══════════════════════════════════════════════════════════════════════════════

xlns16 fp2xlns16(float x) {
    if (x == 0.0f) return SBP_ZERO;
    if (isnan(x))  return SBP_NAN;
    
    int negative = (x < 0.0f);
    float ax = fabsf(x);
    
    // Compute log2(|x|) * 128
    float log2_val = log2f(ax) * SBP_SCALE;
    int16_t log_mag = (int16_t)roundf(log2_val);
    
    // Clamp to representable range
    if (log_mag > 0x3FFF)  log_mag = 0x3FFF;
    if (log_mag < -0x3FFF) log_mag = -0x3FFF;
    
    return sbp_from_sign_log(negative, log_mag);
}

float xlns162fp(xlns16 x) {
    if (sbp_is_zero(x)) return 0.0f;
    if (sbp_is_nan(x))  return nanf("");
    
    int16_t log_mag = sbp_log_mag(x);
    float log2_val = (float)log_mag / SBP_SCALE;
    float magnitude = exp2f(log2_val);
    
    return sbp_is_negative(x) ? -magnitude : magnitude;
}

// ══════════════════════════════════════════════════════════════════════════════
// COMPARISON AND UTILITY
// ══════════════════════════════════════════════════════════════════════════════

inline int xlns16_is_zero(xlns16 x)     { return sbp_is_zero(x); }
inline int xlns16_is_negative(xlns16 x) { return sbp_is_negative(x); }
inline int xlns16_is_positive(xlns16 x) { return sbp_is_positive(x); }

// Comparison using SBP ordering
inline int xlns16_gt(xlns16 a, xlns16 b) {
    int a_neg = sbp_is_negative(a);
    int b_neg = sbp_is_negative(b);
    if (a_neg != b_neg) return b_neg;
    if (sbp_is_zero(a)) return sbp_is_negative(b);
    if (sbp_is_zero(b)) return sbp_is_positive(a);
    int16_t la = sbp_log_mag(a);
    int16_t lb = sbp_log_mag(b);
    if (a_neg) return la < lb;
    return la > lb;
}

inline int xlns16_lt(xlns16 a, xlns16 b) { return xlns16_gt(b, a); }
inline int xlns16_eq(xlns16 a, xlns16 b) { return a == b; }
inline int xlns16_ge(xlns16 a, xlns16 b) { return !xlns16_lt(a, b); }
inline int xlns16_le(xlns16 a, xlns16 b) { return !xlns16_gt(a, b); }
inline xlns16 xlns16_max(xlns16 a, xlns16 b) { return xlns16_gt(a, b) ? a : b; }
inline xlns16 xlns16_min(xlns16 a, xlns16 b) { return xlns16_lt(a, b) ? a : b; }

// ══════════════════════════════════════════════════════════════════════════════
// ADDITIONAL MATH FUNCTIONS
// ══════════════════════════════════════════════════════════════════════════════

inline xlns16 xlns16_square(xlns16 x) { return xlns16_mul(x, x); }

inline xlns16 xlns16_sqrt(xlns16 x) {
    if (sbp_is_zero(x) || sbp_is_negative(x)) return SBP_ZERO;
    int16_t log_mag = sbp_log_mag(x);
    return sbp_from_sign_log(0, log_mag / 2);
}

inline xlns16 xlns16_recip(xlns16 x) {
    if (sbp_is_zero(x)) return SBP_NAN;
    int16_t log_mag = sbp_log_mag(x);
    return sbp_from_sign_log(sbp_is_negative(x), -log_mag);
}

inline xlns16 xlns16_sum(xlns16* x, size_t n) {
    if (n == 0) return SBP_ZERO;
    xlns16 s = x[0];
    for (size_t i = 1; i < n; i++) s = xlns16_add(s, x[i]);
    return s;
}

inline xlns16 xlns16_dot(xlns16* a, xlns16* b, size_t n) {
    if (n == 0) return SBP_ZERO;
    xlns16 s = xlns16_mul(a[0], b[0]);
    for (size_t i = 1; i < n; i++) s = xlns16_add(s, xlns16_mul(a[i], b[i]));
    return s;
}

inline void fp2xlns16_batch(float* src, xlns16* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = fp2xlns16(src[i]);
}

inline void xlns162fp_batch(xlns16* src, float* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = xlns162fp(src[i]);
}

// ══════════════════════════════════════════════════════════════════════════════
// xlns16_float CLASS (C++ wrapper)
// ══════════════════════════════════════════════════════════════════════════════

#ifdef __cplusplus
#include <iostream>

struct xlns16_float {
    xlns16 x;
    
    xlns16_float() : x(SBP_ZERO) {}
    xlns16_float(xlns16 v) : x(v) {}
    
    xlns16_float& operator=(float f) { x = fp2xlns16(f); return *this; }
    
    friend xlns16_float operator+(xlns16_float a, xlns16_float b) { return xlns16_add(a.x, b.x); }
    friend xlns16_float operator-(xlns16_float a, xlns16_float b) { return xlns16_sub(a.x, b.x); }
    friend xlns16_float operator*(xlns16_float a, xlns16_float b) { return xlns16_mul(a.x, b.x); }
    friend xlns16_float operator/(xlns16_float a, xlns16_float b) { return xlns16_div(a.x, b.x); }
    friend xlns16_float operator-(xlns16_float a) { return sbp_neg(a.x); }
    
    friend bool operator==(xlns16_float a, xlns16_float b) { return a.x == b.x; }
    friend bool operator!=(xlns16_float a, xlns16_float b) { return a.x != b.x; }
    friend bool operator<(xlns16_float a, xlns16_float b)  { return xlns16_lt(a.x, b.x); }
    friend bool operator>(xlns16_float a, xlns16_float b)  { return xlns16_gt(a.x, b.x); }
    friend bool operator<=(xlns16_float a, xlns16_float b) { return xlns16_le(a.x, b.x); }
    friend bool operator>=(xlns16_float a, xlns16_float b) { return xlns16_ge(a.x, b.x); }
    
    friend float xlns16_2float(xlns16_float y) { return xlns162fp(y.x); }
    friend xlns16_float float2xlns16_(float f) { xlns16_float r; r.x = fp2xlns16(f); return r; }
    
    friend std::ostream& operator<<(std::ostream& os, xlns16_float y) { return os << xlns162fp(y.x); }
};

inline xlns16_float abs(xlns16_float x)  { return sbp_abs(x.x); }
inline xlns16_float sqrt(xlns16_float x) { return xlns16_sqrt(x.x); }

#endif // __cplusplus
