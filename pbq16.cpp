// pbq16.cpp — Virtual Mirror 16-bit LNS encoding
//
// Index structure:
//
//   code    meaning
//   0x0000  Zero
//   0x0001  largest negative   (dist = 0x7FFF)
//   ...
//   0x7FFE  second smallest negative (dist = 2)
//   0x7FFF  smallest negative  (dist = 1)
//   ──────  axis at 0x7FFF.5
//   0x8000  smallest positive  (dist = 1)   = ~0x7FFF
//   0x8001  second smallest positive (dist = 2)
//   ...
//   0xFFFE  largest positive   (dist = 0x7FFF)  = ~0x0001
//   0xFFFF  NaN                               = ~0x0000
//
// Distance from axis:
//   positive:  dist = code - 0x7FFF          (0x8000 → 1, 0xFFFE → 0x7FFF)
//   negative:  dist = 0x8000 - code          (0x7FFF → 1, 0x0001 → 0x7FFF)
//
// Negation = bitwise NOT (single gate):
//   ~0x7FFF = 0x8000  (smallest neg ↔ smallest pos, dist=1)
//   ~0x0001 = 0xFFFE  (largest neg  ↔ largest pos,  dist=0x7FFF)
//   ~0x0000 = 0xFFFF  (Zero ↔ NaN)
//
// Log magnitude encoding:
//   value = 2^(dist / PBQ_SCALE)
//   PBQ_SCALE = 128  (7 fractional bits)
//
// Arithmetic uses CF-based sb/db (same engine as pbf16_sbp.cpp).

#pragma once
#include <stdint.h>
#include <math.h>

typedef uint16_t pbq16;

// ── Constants ─────────────────────────────────────────────────────────────────

#define PBQ_ZERO      ((pbq16)0x0000u)
#define PBQ_NAN       ((pbq16)0xFFFFu)
#define PBQ_SCALE     128
#define PBQ_MAX_DIST  ((uint16_t)0x7FFFu)

// ── Negation ──────────────────────────────────────────────────────────────────

static inline pbq16 pbq_neg(pbq16 x) { return (pbq16)(~x); }

// ── Predicates ────────────────────────────────────────────────────────────────

static inline int pbq_is_zero(pbq16 x) { return x == PBQ_ZERO; }
static inline int pbq_is_nan(pbq16 x)  { return x == PBQ_NAN; }
static inline int pbq_is_pos(pbq16 x)  { return x >= (pbq16)0x8000u && x != PBQ_NAN; }
static inline int pbq_is_neg(pbq16 x)  { return x > PBQ_ZERO && x < (pbq16)0x8000u; }

// ── Distance from axis ────────────────────────────────────────────────────────

static inline uint16_t pbq_dist(pbq16 x) {
    if (x >= (uint16_t)0x8000u) return (uint16_t)(x - (uint16_t)0x7FFFu);
    return (uint16_t)((uint16_t)0x8000u - x);
}

// ── Construct code from sign + distance ──────────────────────────────────────

static inline pbq16 pbq_make(int negative, uint16_t dist) {
    if (dist == 0)           return PBQ_ZERO;
    if (dist > PBQ_MAX_DIST)  dist = PBQ_MAX_DIST;
    if (negative)            return (pbq16)((uint16_t)0x8000u - dist);
    return (pbq16)((uint16_t)0x7FFFu + dist);
}

// ── Absolute value ────────────────────────────────────────────────────────────

static inline pbq16 pbq_abs(pbq16 x) {
    return pbq_is_neg(x) ? pbq_neg(x) : x;
}

// ── float ↔ pbq16 ──────────────────────────────────────────────────────────────
//
// dist encodes log2(|x|) * PBQ_SCALE, offset so dist=1 maps to the smallest
// representable magnitude.
//
// log2(|x|) ranges over [-PBQ_MAX_DIST/128, +PBQ_MAX_DIST/128].
// We store: dist = (int32_t)(log2(|x|) * PBQ_SCALE) + PBQ_OFFSET
// where PBQ_OFFSET = PBQ_MAX_DIST/2 centres the range.
// dist=1 → log2(|x|) = (1 - PBQ_OFFSET)/128  (smallest magnitude)
// dist=PBQ_MAX_DIST → log2(|x|) = (PBQ_MAX_DIST - PBQ_OFFSET)/128 (largest)

#define PBQ_LOG_OFFSET  ((int32_t)(PBQ_MAX_DIST / 2))   // = 0x3FFF

static pbq16 f2pbq(float x) {
    if (x == 0.0f) return PBQ_ZERO;
    if (x != x)    return PBQ_NAN;
    int neg = (x < 0.0f);
    float ax = neg ? -x : x;
    int32_t dist = (int32_t)(log2f(ax) * (float)PBQ_SCALE + 0.5f) + PBQ_LOG_OFFSET;
    if (dist < 1)           dist = 1;
    if (dist > PBQ_MAX_DIST) dist = PBQ_MAX_DIST;
    return pbq_make(neg, (uint16_t)dist);
}

static float pbq2f(pbq16 x) {
    if (pbq_is_zero(x)) return 0.0f;
    if (pbq_is_nan(x))  return __builtin_nanf("");
    float log2_mag = ((float)(int32_t)pbq_dist(x) - (float)PBQ_LOG_OFFSET) / (float)PBQ_SCALE;
    float mag = exp2f(log2_mag);
    return pbq_is_neg(x) ? -mag : mag;
}

// ── Multiply / divide (exact in log space) ────────────────────────────────────
// log2(a*b) = log2(a) + log2(b)
// dist_a = log2(a)*128 + OFFSET  →  log2(a)*128 = dist_a - OFFSET
// dist_result = (dist_a - OFFSET) + (dist_b - OFFSET) + OFFSET
//             = dist_a + dist_b - OFFSET

static pbq16 pbq_mul(pbq16 a, pbq16 b) {
    if (pbq_is_zero(a) || pbq_is_zero(b)) return PBQ_ZERO;
    if (pbq_is_nan(a)  || pbq_is_nan(b))  return PBQ_NAN;
    int neg = pbq_is_neg(a) ^ pbq_is_neg(b);
    int32_t d = (int32_t)pbq_dist(a) + (int32_t)pbq_dist(b) - PBQ_LOG_OFFSET;
    if (d < 1)           d = 1;
    if (d > PBQ_MAX_DIST) d = PBQ_MAX_DIST;
    return pbq_make(neg, (uint16_t)d);
}

static pbq16 pbq_div(pbq16 a, pbq16 b) {
    if (pbq_is_zero(b) || pbq_is_nan(b)) return PBQ_NAN;
    if (pbq_is_zero(a) || pbq_is_nan(a)) return pbq_is_nan(a) ? PBQ_NAN : PBQ_ZERO;
    int neg = pbq_is_neg(a) ^ pbq_is_neg(b);
    int32_t d = (int32_t)pbq_dist(a) - (int32_t)pbq_dist(b) + PBQ_LOG_OFFSET;
    if (d < 1)           d = 1;
    if (d > PBQ_MAX_DIST) d = PBQ_MAX_DIST;
    return pbq_make(neg, (uint16_t)d);
}

static pbq16 pbq_sq(pbq16 x) { return pbq_mul(x, x); }

static pbq16 pbq_recip(pbq16 x) {
    // log2(1/x) = -log2(x)
    // dist_result = -log2(x)*128 + OFFSET = -(dist_x - OFFSET) + OFFSET = 2*OFFSET - dist_x
    if (pbq_is_zero(x) || pbq_is_nan(x)) return PBQ_NAN;
    int neg = pbq_is_neg(x);
    int32_t d = 2 * PBQ_LOG_OFFSET - (int32_t)pbq_dist(x);
    if (d < 1)           d = 1;
    if (d > PBQ_MAX_DIST) d = PBQ_MAX_DIST;
    return pbq_make(neg, (uint16_t)d);
}

// ── Add / subtract (via float for now; CF sb/db to follow) ───────────────────

static pbq16 pbq_add(pbq16 a, pbq16 b) {
    if (pbq_is_zero(a)) return b;
    if (pbq_is_zero(b)) return a;
    if (pbq_is_nan(a) || pbq_is_nan(b)) return PBQ_NAN;
    return f2pbq(pbq2f(a) + pbq2f(b));
}

static pbq16 pbq_sub(pbq16 a, pbq16 b) { return pbq_add(a, pbq_neg(b)); }
