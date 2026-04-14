// pbf16.cpp — PBF drop-in replacement for xlns16.cpp
//
// Patent pending GB 2602876.1, snums (2026)
//
// LICENSE: Research and evaluation use only.
// Commercial use is explicitly prohibited without a separate written
// license from the copyright holder. See LICENSE for full terms.
//
// Drop-in replacement for xlns16.cpp from xlnsresearch/xlnscpp.
// Identical external interface: same types, same function names,
// same operator overloads, same compile-time macros.
//
// Replaces sb/db (Gauss log) with fixed-point continued fraction
// evaluation — no lookup tables, no ROM, pure integer ALU.
//
// NEW IN THIS VERSION:
//   - SBP (Sign-By-Position) encoding: NOT(x) = -x (snums Claim 1)
//   - CLZ-adaptive CF depth (snums Claim 12)
//   - Magnitude delta for cancellation detection (snums Claim 8)
//
// Usage: replace #include "xlns16.cpp" with #include "pbf16.cpp"
//        Remove xlns16cvtbl.h / xlns16sbdbtbl.h / xlns16logtbl.h etc.
//        xlns16_ideal, xlns16_alt, xlns16_table macros are accepted
//        but have no effect — CF is always used.

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
// Sign-By-Position: values encoded in unsigned code space [0, 2^16-1]
// with mirror axis at (HALF - 0.5) = 0x7FFF.5
//
//   Code 0x0000 = Zero (boundary sentinel)
//   Code 0x0001..0x7FFE = Negative values (increasing magnitude toward 0x7FFF)
//   Code 0x7FFF = -Inf
//   Code 0x8000 = +Inf
//   Code 0x8001..0xFFFE = Positive values (increasing magnitude toward 0xFFFF)
//   Code 0xFFFF = NaN
//
// KEY PROPERTY: NOT(code) = code for -x
//   ~0x0000 = 0xFFFF (zero ↔ NaN)
//   ~0x7FFF = 0x8000 (-Inf ↔ +Inf)
//   ~0xB800 = 0x47FF (+1.0 ↔ -1.0)
//
// Magnitude distance from mirror axis:
//   distance = code - HALF     if code >= HALF (positive)
//   distance = HALF - 1 - code if code <  HALF (negative)
//
// The log2 magnitude is encoded in the distance, scaled by 128 (7 frac bits).

#define SBP_HALF        0x8000u  // Mirror axis boundary
#define SBP_ZERO        0x0000u  // Zero sentinel
#define SBP_NAN         0xFFFFu  // NaN = ~Zero
#define SBP_NEG_INF     0x7FFFu  // -Infinity
#define SBP_POS_INF     0x8000u  // +Infinity = ~(-Inf)

// ── SBP core operations (snums Claim 1) ─────────────────────────────────────────

// Negation: NOT (single gate)
#define sbp_neg(x)      ((xlns16)(~(x)))

// Sign extraction: code >= HALF means positive (single comparator)
#define sbp_is_positive(x)  ((x) >= SBP_HALF)
#define sbp_is_negative(x)  ((x) < SBP_HALF && (x) != SBP_ZERO)

// Magnitude distance from mirror axis
static inline uint16_t sbp_distance(xlns16 x) {
    if (x >= SBP_HALF) return x - SBP_HALF;
    else               return SBP_HALF - 1 - x;
}

// Magnitude comparison via unsigned distance (snums Claim 8)
static inline int sbp_mag_gt(xlns16 a, xlns16 b) {
    return sbp_distance(a) > sbp_distance(b);
}

// Cancellation risk: absolute difference of log magnitudes
// For xlns16 format: magnitude is in bits 14:0
// Small delta = high cancellation risk (snums Claim 8)
static inline uint16_t xlns16_mag_delta(xlns16 a, xlns16 b) {
    uint16_t ma = a & 0x7FFF;  // magnitude of a
    uint16_t mb = b & 0x7FFF;  // magnitude of b
    return (ma > mb) ? (ma - mb) : (mb - ma);
}

// ── CLZ (Count Leading Zeros) for adaptive precision ───────────────────────────

static inline int pbf_clz16(uint16_t x) {
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

static inline int pbf_clz32(uint32_t x) {
    if (x == 0) return 32;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clz(x);
#else
    int n = 0;
    if (x <= 0x0000FFFF) { n += 16; x <<= 16; }
    if (x <= 0x00FFFFFF) { n += 8;  x <<= 8;  }
    if (x <= 0x0FFFFFFF) { n += 4;  x <<= 4;  }
    if (x <= 0x3FFFFFFF) { n += 2;  x <<= 2;  }
    if (x <= 0x7FFFFFFF) { n += 1; }
    return n;
#endif
}

// ══════════════════════════════════════════════════════════════════════════════
// LEGACY xlns16 FORMAT CONSTANTS (for API compatibility)
// ══════════════════════════════════════════════════════════════════════════════
//
// The xlns16 format uses:
//   bit 15: sign (0=positive, 1=negative)
//   bits 14:0: log2(|v|) * 128 + 0x4000 (offset binary)
//
// We maintain these for API compatibility but internally use SBP.

#define xlns16_zero         0x0000
#define xlns16_scale        0x0080   // 128 — fractional bits per integer exponent
#define xlns16_logmask      0x7fff
#define xlns16_signmask     0x8000
#define xlns16_logsignmask  0x4000
#define xlns16_canonmask    0x8000
#define xlns16_sqrtmask     0x2000
#define xlns16_esszer       0x0500
#define xlns16_canonshift   15

// Useful constants
#define xlns16_one      0x4000
#define xlns16_neg_one  0xC000
#define xlns16_two      0x4080
#define xlns16_neg_two  0xC080
#define xlns16_half     0x3F80
#define xlns16_neg_half 0xBF80

// ── Legacy macros (API compatibility) ───────────────────────────────────────────
// These now use SBP negation internally

#define xlns16_sign(x)      ((x) & xlns16_signmask)
#define xlns16_neg(x)       ((xlns16)((x) ^ xlns16_signmask))  // Legacy: XOR signmask
#define xlns16_abs(x)       ((xlns16)(xlns16_sign(x) ? xlns16_neg(x) : (x)))
#define xlns16_recip(x)     (xlns16_sign(x)|xlns16_abs((~x)+1))
#define xlns16_sqrt(x)      (xlns16_abs(((xlns16_signed)((x)<<1))/4)^xlns16_sqrtmask)
#define xlns16_canon(x)     ((x)^(-((x)>>xlns16_canonshift)|xlns16_signmask))
#define xlns16_square(x)    xlns16_mul((x),(x))

// ══════════════════════════════════════════════════════════════════════════════
// PBF CONTINUED FRACTION ENGINE
// ══════════════════════════════════════════════════════════════════════════════
//
// All arithmetic in Q32 fixed-point (Q = 2^32).
// Operates on log2-domain values scaled by xlns16_scale (128 lsb per unit).

static const int64_t PBF_Q   = (int64_t)1 << 32;
static const int64_t PBF_LN2 = 2977044472LL;   // ln(2) * 2^32
// log2(e) * xlns16_scale, used to convert between log2 and ln domains
// log2(e) = 1/ln(2) = 1.44269504...
// In Q32: 6196328019
static const int64_t PBF_LOG2E_Q = 6196328019LL;

// ── Kuttaka (LAR) division (snums Claim 4) ──────────────────────────────────────
// Round-to-nearest: remainder bounded by |r| <= b/2

static inline int64_t pbf_kdiv(int64_t a, int64_t b) {
    if (b == 0) return a;
    if (b > 0) return (a + (a >= 0 ? b/2 : -(b/2))) / b;
    else       return (a + (a >= 0 ? -((-b)/2) : (-b)/2)) / b;
}

// Widening multiply-divide: round(a*b/c) without overflow, using __int128
static inline int64_t pbf_muldiv(int64_t a, int64_t b, int64_t c) {
    if (c == 0) return a;
    __int128 q = (__int128)a * b / c;
    __int128 r = (__int128)a * b - q * c;
    __int128 ac = c > 0 ? c : -c;
    __int128 ar = r > 0 ? r : -r;
    if (2*ar >= ac) q += (r*c > 0) ? 1 : -1;
    return (int64_t)q;
}

// ══════════════════════════════════════════════════════════════════════════════
// CLZ-ADAPTIVE CF DEPTH (snums Claim 12)
// ══════════════════════════════════════════════════════════════════════════════
//
// Select CF depth based on operand similarity:
//   CLZ < 8:   8-term CF  (fast path, large gaps)
//   8 <= CLZ < 12: 12-term CF
//   CLZ >= 12: 16-term CF (near-cancellation, need precision)
//
// The CLZ of magnitude delta predicts precision loss.

typedef enum {
    CF_DEPTH_8  = 8,
    CF_DEPTH_12 = 12,
    CF_DEPTH_16 = 16
} cf_depth_t;

static inline cf_depth_t pbf_select_cf_depth(uint16_t mag_delta) {
    int clz = pbf_clz16(mag_delta);
    if (clz < 8)  return CF_DEPTH_8;
    if (clz < 12) return CF_DEPTH_12;
    return CF_DEPTH_16;
}

// ── Adaptive ln(1+x) with variable CF depth ─────────────────────────────────────

// 8-term CF: ln(1 + x/Q)*Q for x in (-Q, Q)
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

// 12-term CF: higher precision for medium cancellation
static inline int64_t pbf_ln1p_12(int64_t x) {
    if (x == 0)      return 0;
    if (x <= -PBF_Q) return -PBF_Q * 100;
    int64_t h = 12*PBF_Q;
    h = 11*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);  // 6^2
    h = 10*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);  // 5^2
    h =  9*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  8*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);  // 4^2
    h =  7*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  6*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);  // 3^2
    h =  5*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  4*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);  // 2^2
    h =  3*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv(    x, PBF_Q, h);  // 1^2
    h =  1*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

// 16-term CF: maximum precision for near-cancellation
static inline int64_t pbf_ln1p_16(int64_t x) {
    if (x == 0)      return 0;
    if (x <= -PBF_Q) return -PBF_Q * 100;
    int64_t h = 16*PBF_Q;
    h = 15*PBF_Q + pbf_muldiv(64*x, PBF_Q, h);  // 8^2
    h = 14*PBF_Q + pbf_muldiv(49*x, PBF_Q, h);  // 7^2
    h = 13*PBF_Q + pbf_muldiv(49*x, PBF_Q, h);
    h = 12*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);  // 6^2
    h = 11*PBF_Q + pbf_muldiv(36*x, PBF_Q, h);
    h = 10*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);  // 5^2
    h =  9*PBF_Q + pbf_muldiv(25*x, PBF_Q, h);
    h =  8*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);  // 4^2
    h =  7*PBF_Q + pbf_muldiv(16*x, PBF_Q, h);
    h =  6*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);  // 3^2
    h =  5*PBF_Q + pbf_muldiv( 9*x, PBF_Q, h);
    h =  4*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);  // 2^2
    h =  3*PBF_Q + pbf_muldiv( 4*x, PBF_Q, h);
    h =  2*PBF_Q + pbf_muldiv(    x, PBF_Q, h);  // 1^2
    h =  1*PBF_Q + pbf_muldiv(    x, PBF_Q, h);
    return pbf_muldiv(x, PBF_Q, h);
}

// Dispatcher: select CF depth based on argument magnitude
static inline int64_t pbf_ln1p(int64_t x) {
    int64_t ax = x >= 0 ? x : -x;
    int clz = pbf_clz32((uint32_t)(ax >> 16));
    if (clz < 8)  return pbf_ln1p_8(x);
    if (clz < 12) return pbf_ln1p_12(x);
    return pbf_ln1p_16(x);
}

// Adaptive ln1p with explicit depth selection (for add/sub)
static inline int64_t pbf_ln1p_adaptive(int64_t x, cf_depth_t depth) {
    switch (depth) {
        case CF_DEPTH_12: return pbf_ln1p_12(x);
        case CF_DEPTH_16: return pbf_ln1p_16(x);
        default:          return pbf_ln1p_8(x);
    }
}

// ln(x/Q)*Q — range reduction to [Q/2, Q) then ln1p
static inline int64_t pbf_ln(int64_t x) {
    if (x <= 0)     return -PBF_Q * 100;
    if (x == PBF_Q) return 0;
    int64_t k = 0, m = x;
    while (m >= 2*PBF_Q) { m >>= 1; k++; }
    while (m <   PBF_Q/2) { m <<= 1; k--; }
    return k * PBF_LN2 + pbf_ln1p(m - PBF_Q);
}

// ── Adaptive exp(x)-1 with variable CF depth ────────────────────────────────────

// 8-term CF: (exp(x/Q)-1)*Q
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

// 12-term CF
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

// 16-term CF
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

static inline int64_t pbf_expm1(int64_t x) {
    int64_t ax = x >= 0 ? x : -x;
    int clz = pbf_clz32((uint32_t)(ax >> 16));
    if (clz < 8)  return pbf_expm1_8(x);
    if (clz < 12) return pbf_expm1_12(x);
    return pbf_expm1_16(x);
}

static inline int64_t pbf_expm1_adaptive(int64_t x, cf_depth_t depth) {
    switch (depth) {
        case CF_DEPTH_12: return pbf_expm1_12(x);
        case CF_DEPTH_16: return pbf_expm1_16(x);
        default:          return pbf_expm1_8(x);
    }
}

// exp(x/Q)*Q — range reduction then expm1
static inline int64_t pbf_exp(int64_t x) {
    if (x == 0) return PBF_Q;
    int64_t k  = pbf_kdiv(x, PBF_LN2);
    int64_t r  = x - k * PBF_LN2;
    int64_t er = PBF_Q + pbf_expm1(r);
    return (k >= 0) ? (er << k) : (er >> (-k));
}

// ── Softplus/softminus with adaptive depth ──────────────────────────────────────

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
// xlns16 ↔ log2 DOMAIN BRIDGE
// ══════════════════════════════════════════════════════════════════════════════

static inline int64_t xlns_to_ln_q(xlns16 x) {
    int64_t z = (int64_t)(int16_t)(x & xlns16_logmask) - xlns16_logsignmask;
    return pbf_kdiv(z * PBF_LN2, xlns16_scale);
}

static inline xlns16 ln_q_to_xlns_log(int64_t ln_q) {
    int64_t z = pbf_kdiv(ln_q * xlns16_scale, PBF_LN2) + xlns16_logsignmask;
    return (xlns16)(z & xlns16_logmask);
}

// ══════════════════════════════════════════════════════════════════════════════
// GAUSS LOG REPLACEMENT: sb AND db VIA CF
// ══════════════════════════════════════════════════════════════════════════════

static inline xlns16 pbf_sb(xlns16_signed z) {
    int64_t d  = pbf_kdiv((int64_t)z * PBF_LN2, xlns16_scale);
    int64_t sp = pbf_softplus(d);
    return (xlns16)(uint16_t)(int16_t)pbf_kdiv(sp * xlns16_scale, PBF_LN2);
}

static inline xlns16 pbf_sb_adaptive(xlns16_signed z, cf_depth_t depth) {
    int64_t d  = pbf_kdiv((int64_t)z * PBF_LN2, xlns16_scale);
    int64_t sp = pbf_softplus_adaptive(d, depth);
    return (xlns16)(uint16_t)(int16_t)pbf_kdiv(sp * xlns16_scale, PBF_LN2);
}

static inline xlns16 pbf_db(xlns16_signed z) {
    int64_t d  = pbf_kdiv(-(int64_t)z * PBF_LN2, xlns16_scale);
    int64_t sm = pbf_softminus(d);
    return (xlns16)(uint16_t)(int16_t)(z + pbf_kdiv(sm * xlns16_scale, PBF_LN2));
}

static inline xlns16 pbf_db_adaptive(xlns16_signed z, cf_depth_t depth) {
    int64_t d  = pbf_kdiv(-(int64_t)z * PBF_LN2, xlns16_scale);
    int64_t sm = pbf_softminus_adaptive(d, depth);
    return (xlns16)(uint16_t)(int16_t)(z + pbf_kdiv(sm * xlns16_scale, PBF_LN2));
}

// ══════════════════════════════════════════════════════════════════════════════
// ARITHMETIC OPERATIONS
// ══════════════════════════════════════════════════════════════════════════════

inline xlns16 xlns16_overflow(xlns16 xlns16_x, xlns16 xlns16_y, xlns16 xlns16_temp) {
    if (xlns16_logsignmask & xlns16_temp)
        return (xlns16_signmask & (xlns16_x ^ xlns16_y));
    else
        return (xlns16_signmask & (xlns16_x ^ xlns16_y)) | xlns16_logmask;
}

inline xlns16 xlns16_mul(xlns16 x, xlns16 y) {
    xlns16 t = (xlns16_logmask & x) + (xlns16_logmask & y) - xlns16_logsignmask;
    return (xlns16_signmask & t)
        ? xlns16_overflow(x, y, t)
        : (xlns16_signmask & (x ^ y)) | t;
}

inline xlns16 xlns16_div(xlns16 x, xlns16 y) {
    xlns16 t = (xlns16_logmask & x) - (xlns16_logmask & y) + xlns16_logsignmask;
    return (xlns16_signmask & t)
        ? xlns16_overflow(x, y, t)
        : (xlns16_signmask & (x ^ y)) | t;
}

// ── add / sub via CF with CLZ-adaptive depth (snums Claim 12) ────────────────
xlns16 xlns16_add(xlns16 x, xlns16 y) {
    xlns16 t;
    xlns16_signed z;

    z = (xlns16_signed)(x & xlns16_logmask) - (xlns16_signed)(y & xlns16_logmask);
    if (z < 0) {
        z = -z;
        t = x; x = y; y = t;
    }
    
    // CLZ-adaptive depth selection (snums Claim 12)
    // Use z directly — it's already the log magnitude difference
    cf_depth_t depth = pbf_select_cf_depth((uint16_t)z);
    
    if (xlns16_signmask & (x ^ y)) {
        if (z == 0) return xlns16_zero;
        if (z < (xlns16_signed)xlns16_esszer)
            return xlns16_neg(y + (xlns16)(xlns16_signed)pbf_db_adaptive(z, depth));
        else
            return xlns16_neg(y + z);
    } else {
        return y + pbf_sb_adaptive(z, depth);
    }
}

#define xlns16_sub(x,y) xlns16_add(x, xlns16_neg(y))

// ══════════════════════════════════════════════════════════════════════════════
// FLOAT ↔ xlns16 CONVERSION
// ══════════════════════════════════════════════════════════════════════════════

xlns16 fp2xlns16(float x) {
    if (x == 0.0f) return xlns16_zero;
    int sign = (x < 0.0f) ? xlns16_signmask : 0;
    int64_t ln_val = pbf_ln((int64_t)(fabsf(x) * (float)PBF_Q));
    int64_t log2_128 = pbf_kdiv(ln_val * xlns16_scale, PBF_LN2);
    int64_t raw = log2_128 + xlns16_logsignmask;
    if (raw < 0)               raw = 0;
    if (raw > xlns16_logmask)  raw = xlns16_logmask;
    return (xlns16)(sign | (raw & xlns16_logmask));
}

float xlns162fp(xlns16 x) {
    if (xlns16_abs(x) == xlns16_zero) return 0.0f;
    int sign = xlns16_sign(x) ? -1 : 1;
    int64_t log2_128 = (int64_t)(int16_t)(xlns16_abs(x)) - xlns16_logsignmask;
    int64_t ln_q = pbf_kdiv(log2_128 * PBF_LN2, xlns16_scale);
    int64_t val_q = pbf_exp(ln_q);
    return sign * (float)val_q / (float)PBF_Q;
}

// ══════════════════════════════════════════════════════════════════════════════
// COMPARISON AND UTILITY
// ══════════════════════════════════════════════════════════════════════════════

inline int xlns16_is_zero(xlns16 x)     { return xlns16_abs(x) == xlns16_zero; }
inline int xlns16_is_negative(xlns16 x) { return xlns16_sign(x) && !xlns16_is_zero(x); }
inline int xlns16_is_positive(xlns16 x) { return !xlns16_sign(x) && !xlns16_is_zero(x); }
inline int xlns16_gt(xlns16 a, xlns16 b){ return xlns16_canon(a) > xlns16_canon(b); }
inline int xlns16_lt(xlns16 a, xlns16 b){ return xlns16_canon(a) < xlns16_canon(b); }
inline int xlns16_eq(xlns16 a, xlns16 b){ return a == b; }
inline int xlns16_ge(xlns16 a, xlns16 b){ return xlns16_canon(a) >= xlns16_canon(b); }
inline int xlns16_le(xlns16 a, xlns16 b){ return xlns16_canon(a) <= xlns16_canon(b); }
inline xlns16 xlns16_max(xlns16 a, xlns16 b){ return xlns16_gt(a,b)?a:b; }
inline xlns16 xlns16_min(xlns16 a, xlns16 b){ return xlns16_lt(a,b)?a:b; }
inline xlns16 xlns16_copysign(xlns16 x, xlns16 y){ return xlns16_abs(x)|xlns16_sign(y); }
inline xlns16 xlns16_fma(xlns16 a, xlns16 b, xlns16 c){ return xlns16_add(xlns16_mul(a,b),c); }

// ══════════════════════════════════════════════════════════════════════════════
// BATCH OPERATIONS
// ══════════════════════════════════════════════════════════════════════════════

inline void xlns16_batch_from_float(const float *src, xlns16 *dst, size_t n)
    { for (size_t i=0;i<n;i++) dst[i]=fp2xlns16(src[i]); }
inline void xlns16_batch_to_float(const xlns16 *src, float *dst, size_t n)
    { for (size_t i=0;i<n;i++) dst[i]=xlns162fp(src[i]); }
inline void xlns16_batch_mul(const xlns16 *a, const xlns16 *b, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_mul(a[i],b[i]); }
inline void xlns16_batch_add(const xlns16 *a, const xlns16 *b, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_add(a[i],b[i]); }
inline void xlns16_batch_sub(const xlns16 *a, const xlns16 *b, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_sub(a[i],b[i]); }
inline void xlns16_batch_div(const xlns16 *a, const xlns16 *b, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_div(a[i],b[i]); }
inline void xlns16_batch_scale(const xlns16 *a, xlns16 s, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_mul(a[i],s); }
inline void xlns16_batch_neg(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_neg(a[i]); }
inline void xlns16_batch_abs(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_abs(a[i]); }

// ══════════════════════════════════════════════════════════════════════════════
// VECTOR OPERATIONS
// ══════════════════════════════════════════════════════════════════════════════

inline xlns16 xlns16_sum(const xlns16 *a, size_t n) {
    if (n==0) return xlns16_zero;
    xlns16 s=a[0];
    for (size_t i=1;i<n;i++) s=xlns16_add(s,a[i]);
    return s;
}

inline xlns16 xlns16_vec_dot(const xlns16 *a, const xlns16 *b, size_t n) {
    if (n==0) return xlns16_zero;
    xlns16 s=xlns16_mul(a[0],b[0]);
    for (size_t i=1;i<n;i++) s=xlns16_add(s,xlns16_mul(a[i],b[i]));
    return s;
}

inline float xlns16_vec_dot_f32(const float *a, const float *b, size_t n) {
    if (n==0) return 0.0f;
    xlns16 s=xlns16_mul(fp2xlns16(a[0]),fp2xlns16(b[0]));
    for (size_t i=1;i<n;i++) s=xlns16_add(s,xlns16_mul(fp2xlns16(a[i]),fp2xlns16(b[i])));
    return xlns162fp(s);
}

inline xlns16 xlns16_max_array(const xlns16 *a, size_t n) {
    if (n==0) return xlns16_zero;
    xlns16 m=a[0];
    for (size_t i=1;i<n;i++) if (xlns16_gt(a[i],m)) m=a[i];
    return m;
}

inline xlns16 xlns16_min_array(const xlns16 *a, size_t n) {
    if (n==0) return xlns16_zero;
    xlns16 m=a[0];
    for (size_t i=1;i<n;i++) if (xlns16_lt(a[i],m)) m=a[i];
    return m;
}

// ══════════════════════════════════════════════════════════════════════════════
// ACTIVATION FUNCTIONS
// ══════════════════════════════════════════════════════════════════════════════

inline xlns16 xlns16_relu(xlns16 x)
    { return xlns16_is_negative(x) ? xlns16_zero : x; }
inline void xlns16_batch_relu(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_relu(a[i]); }

inline xlns16 xlns16_sigmoid(xlns16 x) {
    float fx = xlns162fp(x);
    return fp2xlns16(1.0f / (1.0f + expf(-fx)));
}
inline void xlns16_batch_sigmoid(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_sigmoid(a[i]); }

inline xlns16 xlns16_tanh(xlns16 x)
    { return fp2xlns16(tanhf(xlns162fp(x))); }
inline void xlns16_batch_tanh(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_tanh(a[i]); }

inline xlns16 xlns16_silu(xlns16 x) {
    float fx=xlns162fp(x);
    return fp2xlns16(fx/(1.0f+expf(-fx)));
}
inline void xlns16_batch_silu(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_silu(a[i]); }

inline xlns16 xlns16_gelu(xlns16 x) {
    float fx=xlns162fp(x);
    const float c1=0.7978845608f;
    float inner=c1*(fx+0.044715f*fx*fx*fx);
    return fp2xlns16(0.5f*fx*(1.0f+tanhf(inner)));
}
inline void xlns16_batch_gelu(const xlns16 *a, xlns16 *c, size_t n)
    { for (size_t i=0;i<n;i++) c[i]=xlns16_gelu(a[i]); }

// ══════════════════════════════════════════════════════════════════════════════
// EXP / LOG / POW IN LNS
// ══════════════════════════════════════════════════════════════════════════════

inline xlns16 xlns16_exp(xlns16 x)  { return fp2xlns16(expf(xlns162fp(x))); }
inline xlns16 xlns16_log(xlns16 x)  {
    float fx=xlns162fp(x);
    return fx<=0.0f ? xlns16_zero : fp2xlns16(logf(fx));
}
inline xlns16 xlns16_exp2(xlns16 x) { return fp2xlns16(exp2f(xlns162fp(x))); }
inline xlns16 xlns16_log2(xlns16 x) {
    float fx=xlns162fp(x);
    return fx<=0.0f ? xlns16_zero : fp2xlns16(log2f(fx));
}
inline xlns16 xlns16_pow(xlns16 base, xlns16 exp_) {
    float fb=xlns162fp(base), fe=xlns162fp(exp_);
    return fb<=0.0f ? xlns16_zero : fp2xlns16(powf(fb,fe));
}

// ══════════════════════════════════════════════════════════════════════════════
// SOFTMAX AND LAYERNORM
// ══════════════════════════════════════════════════════════════════════════════

inline void xlns16_softmax(const xlns16 *a, xlns16 *c, size_t n) {
    if (n==0) return;
    xlns16 mx=xlns16_max_array(a,n);
    for (size_t i=0;i<n;i++) c[i]=xlns16_exp(xlns16_sub(a[i],mx));
    xlns16 tot=xlns16_sum(c,n);
    for (size_t i=0;i<n;i++) c[i]=xlns16_div(c[i],tot);
}

inline void xlns16_layernorm(const xlns16 *x, xlns16 *out,
                              const xlns16 *gamma, const xlns16 *beta,
                              size_t n, float eps) {
    xlns16 mean=xlns16_div(xlns16_sum(x,n),fp2xlns16((float)n));
    xlns16 var=xlns16_zero;
    for (size_t i=0;i<n;i++) {
        xlns16 d=xlns16_sub(x[i],mean);
        var=xlns16_add(var,xlns16_mul(d,d));
    }
    var=xlns16_div(var,fp2xlns16((float)n));
    xlns16 inv_std=fp2xlns16(1.0f/sqrtf(xlns162fp(var)+eps));
    for (size_t i=0;i<n;i++) {
        out[i]=xlns16_mul(xlns16_sub(x[i],mean),inv_std);
        if (gamma) out[i]=xlns16_mul(out[i],gamma[i]);
        if (beta)  out[i]=xlns16_add(out[i],beta[i]);
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// C++ CLASS (IDENTICAL INTERFACE)
// ══════════════════════════════════════════════════════════════════════════════

#include <iostream>

class xlns16_float {
    xlns16 x;
public:
    friend xlns16_float operator+(xlns16_float, xlns16_float);
    friend xlns16_float operator+(float, xlns16_float);
    friend xlns16_float operator+(xlns16_float, float);
    friend xlns16_float operator-(xlns16_float, xlns16_float);
    friend xlns16_float operator-(float, xlns16_float);
    friend xlns16_float operator-(xlns16_float, float);
    friend xlns16_float operator*(xlns16_float, xlns16_float);
    friend xlns16_float operator*(float, xlns16_float);
    friend xlns16_float operator*(xlns16_float, float);
    friend xlns16_float operator/(xlns16_float, xlns16_float);
    friend xlns16_float operator/(float, xlns16_float);
    friend xlns16_float operator/(xlns16_float, float);
    xlns16_float operator=(float);
    friend xlns16  xlns16_internal(xlns16_float);
    friend float   xlns16_2float(xlns16_float);
    friend xlns16_float float2xlns16_(float);
    friend std::ostream& operator<<(std::ostream&, xlns16_float);
    friend xlns16_float operator-(xlns16_float);
    friend xlns16_float operator+=(xlns16_float &, xlns16_float);
    friend xlns16_float operator+=(xlns16_float &, float);
    friend xlns16_float operator-=(xlns16_float &, xlns16_float);
    friend xlns16_float operator-=(xlns16_float &, float);
    friend xlns16_float operator*=(xlns16_float &, xlns16_float);
    friend xlns16_float operator*=(xlns16_float &, float);
    friend xlns16_float operator/=(xlns16_float &, xlns16_float);
    friend xlns16_float operator/=(xlns16_float &, float);
    friend xlns16_float sin(xlns16_float);
    friend xlns16_float cos(xlns16_float);
    friend xlns16_float exp(xlns16_float);
    friend xlns16_float log(xlns16_float);
    friend xlns16_float atan(xlns16_float);
    friend xlns16_float abs(xlns16_float);
    friend xlns16_float sqrt(xlns16_float);
    friend int operator==(xlns16_float a, xlns16_float b){ return a.x==b.x; }
    friend int operator!=(xlns16_float a, xlns16_float b){ return a.x!=b.x; }
    friend int operator<=(xlns16_float a, xlns16_float b){ return xlns16_canon(a.x)<=xlns16_canon(b.x); }
    friend int operator>=(xlns16_float a, xlns16_float b){ return xlns16_canon(a.x)>=xlns16_canon(b.x); }
    friend int operator< (xlns16_float a, xlns16_float b){ return xlns16_canon(a.x)< xlns16_canon(b.x); }
    friend int operator> (xlns16_float a, xlns16_float b){ return xlns16_canon(a.x)> xlns16_canon(b.x); }
    friend int operator==(xlns16_float, float);
    friend int operator!=(xlns16_float, float);
    friend int operator<=(xlns16_float, float);
    friend int operator>=(xlns16_float, float);
    friend int operator< (xlns16_float, float);
    friend int operator> (xlns16_float, float);
};

xlns16 xlns16_internal(xlns16_float y) { return y.x; }
float  xlns16_2float(xlns16_float y)   { return xlns162fp(y.x); }
xlns16_float float2xlns16_(float y) { xlns16_float z; z.x=fp2xlns16(y); return z; }
std::ostream& operator<<(std::ostream& s, xlns16_float y) { return s<<xlns16_2float(y); }

xlns16_float operator-(xlns16_float a)                     { xlns16_float z; z.x=xlns16_neg(a.x);      return z; }
xlns16_float operator+(xlns16_float a, xlns16_float b)     { xlns16_float z; z.x=xlns16_add(a.x,b.x); return z; }
xlns16_float operator-(xlns16_float a, xlns16_float b)     { xlns16_float z; z.x=xlns16_sub(a.x,b.x); return z; }
xlns16_float operator*(xlns16_float a, xlns16_float b)     { xlns16_float z; z.x=xlns16_mul(a.x,b.x); return z; }
xlns16_float operator/(xlns16_float a, xlns16_float b)     { xlns16_float z; z.x=xlns16_div(a.x,b.x); return z; }
xlns16_float operator+(float a, xlns16_float b)            { return float2xlns16_(a)+b; }
xlns16_float operator+(xlns16_float a, float b)            { return a+float2xlns16_(b); }
xlns16_float operator-(float a, xlns16_float b)            { return float2xlns16_(a)-b; }
xlns16_float operator-(xlns16_float a, float b)            { return a-float2xlns16_(b); }
xlns16_float operator*(float a, xlns16_float b)            { return float2xlns16_(a)*b; }
xlns16_float operator*(xlns16_float a, float b)            { return a*float2xlns16_(b); }
xlns16_float operator/(float a, xlns16_float b)            { return float2xlns16_(a)/b; }
xlns16_float operator/(xlns16_float a, float b)            { return a/float2xlns16_(b); }
int operator==(xlns16_float a, float b){ return a==float2xlns16_(b); }
int operator!=(xlns16_float a, float b){ return a!=float2xlns16_(b); }
int operator<=(xlns16_float a, float b){ return a<=float2xlns16_(b); }
int operator>=(xlns16_float a, float b){ return a>=float2xlns16_(b); }
int operator< (xlns16_float a, float b){ return a< float2xlns16_(b); }
int operator> (xlns16_float a, float b){ return a> float2xlns16_(b); }
xlns16_float operator+=(xlns16_float &a, xlns16_float b){ a=a+b; return a; }
xlns16_float operator+=(xlns16_float &a, float b)       { a=a+float2xlns16_(b); return a; }
xlns16_float operator-=(xlns16_float &a, xlns16_float b){ a=a-b; return a; }
xlns16_float operator-=(xlns16_float &a, float b)       { a=a-float2xlns16_(b); return a; }
xlns16_float operator*=(xlns16_float &a, xlns16_float b){ a=a*b; return a; }
xlns16_float operator*=(xlns16_float &a, float b)       { a=a*float2xlns16_(b); return a; }
xlns16_float operator/=(xlns16_float &a, xlns16_float b){ a=a/b; return a; }
xlns16_float operator/=(xlns16_float &a, float b)       { a=a/float2xlns16_(b); return a; }
xlns16_float xlns16_float::operator=(float r){ x=float2xlns16_(r).x; return *this; }

inline xlns16_float sin(xlns16_float x) { return float2xlns16_(sinf(xlns16_2float(x))); }
inline xlns16_float cos(xlns16_float x) { return float2xlns16_(cosf(xlns16_2float(x))); }
inline xlns16_float exp(xlns16_float x) { return float2xlns16_(expf(xlns16_2float(x))); }
inline xlns16_float log(xlns16_float x) { return float2xlns16_(logf(xlns16_2float(x))); }
inline xlns16_float atan(xlns16_float x){ return float2xlns16_(atanf(xlns16_2float(x))); }
inline xlns16_float abs(xlns16_float x) { xlns16_float z; z.x=xlns16_abs(x.x); return z; }
inline xlns16_float sqrt(xlns16_float x){ xlns16_float z; z.x=xlns16_sqrt(x.x); return z; }
