// pbf_xlns.cpp — Drop-in xlns16 API shim backed by PBF16.
//
// Provides the xlns16_* function set (fp2xlns16, xlns16_add, xlns16_mul, ...)
// and constants (xlns16_one, xlns16_zero, ...) so code written against the
// xlns16 API can run against the parameterized PBF16 format.
//
// NOTE: the underlying bit patterns differ from xlns16 (PBF uses
// zero-at-0 / inf-at-max / sign-by-position, not log-with-sign-bit).
// Numerical results match; bit-exact encodings do not.

#pragma once

#include "pbf_batch.cpp"
#include <cmath>

// ── Type aliases matching xlns16 naming ──────────────────────────────────────

typedef uint16_t xlns16;
typedef int16_t  xlns16_signed;

// ── Lazy-initialized PBF16 context (used by every xlns16_* entrypoint) ───────

static inline const pbf_t* pbf_xlns_ctx() {
    static pbf_t s_ctx = []{ pbf_t p; pbf16_init(&p); return p; }();
    return &s_ctx;
}

// Convert between xlns16 (uint16_t) and PBF code (uint32_t).
// For PBF16, max_code = 0xFFFF so codes fit in uint16_t exactly.
static inline xlns16   pbf_to_x16(uint32_t c) { return (xlns16)c; }
static inline uint32_t x16_to_pbf(xlns16 x)   { return (uint32_t)x; }

// ── Scalar conversions ───────────────────────────────────────────────────────

static inline xlns16 fp2xlns16(float x) {
    return pbf_to_x16(pbf_encode(pbf_xlns_ctx(), (double)x));
}

static inline float xlns162fp(xlns16 x) {
    return (float)pbf_decode(pbf_xlns_ctx(), x16_to_pbf(x));
}

// ── Constants (computed at static-init time from PBF16 encoding) ─────────────

static const xlns16 xlns16_zero     = pbf_to_x16(pbf_encode(pbf_xlns_ctx(),  0.0));
static const xlns16 xlns16_one      = pbf_to_x16(pbf_encode(pbf_xlns_ctx(),  1.0));
static const xlns16 xlns16_neg_one  = pbf_to_x16(pbf_encode(pbf_xlns_ctx(), -1.0));
static const xlns16 xlns16_two      = pbf_to_x16(pbf_encode(pbf_xlns_ctx(),  2.0));
static const xlns16 xlns16_neg_two  = pbf_to_x16(pbf_encode(pbf_xlns_ctx(), -2.0));
static const xlns16 xlns16_half     = pbf_to_x16(pbf_encode(pbf_xlns_ctx(),  0.5));
static const xlns16 xlns16_neg_half = pbf_to_x16(pbf_encode(pbf_xlns_ctx(), -0.5));

// ── Core arithmetic ──────────────────────────────────────────────────────────

static inline xlns16 xlns16_neg(xlns16 x) {
    return pbf_to_x16(pbf_neg(pbf_xlns_ctx(), x16_to_pbf(x)));
}

static inline xlns16 xlns16_add(xlns16 a, xlns16 b) {
    return pbf_to_x16(pbf_add(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)));
}

static inline xlns16 xlns16_sub(xlns16 a, xlns16 b) {
    return pbf_to_x16(pbf_sub(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)));
}

static inline xlns16 xlns16_mul(xlns16 a, xlns16 b) {
    return pbf_to_x16(pbf_mul(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)));
}

static inline xlns16 xlns16_div(xlns16 a, xlns16 b) {
    return pbf_to_x16(pbf_div(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)));
}

static inline xlns16 xlns16_abs(xlns16 x) {
    const pbf_t* p = pbf_xlns_ctx();
    uint32_t c = x16_to_pbf(x);
    if (c == 0 || c >= p->mid) return x;
    return pbf_to_x16(pbf_neg(p, c));
}

static inline xlns16 xlns16_recip(xlns16 x) {
    return pbf_to_x16(pbf_div(pbf_xlns_ctx(), x16_to_pbf(xlns16_one), x16_to_pbf(x)));
}

static inline xlns16 xlns16_square(xlns16 x) { return xlns16_mul(x, x); }

static inline xlns16 xlns16_sqrt(xlns16 x) {
    float f = xlns162fp(x);
    if (f <= 0.0f) return xlns16_zero;
    return fp2xlns16(sqrtf(f));
}

static inline xlns16 xlns16_fma(xlns16 a, xlns16 b, xlns16 c) {
    return xlns16_add(xlns16_mul(a, b), c);
}

// ── Predicates & sign ────────────────────────────────────────────────────────

static inline int xlns16_is_zero(xlns16 x) { return x == xlns16_zero; }

static inline int xlns16_is_negative(xlns16 x) {
    const pbf_t* p = pbf_xlns_ctx();
    uint32_t c = x16_to_pbf(x);
    return (c != 0) && (c < p->mid);
}

static inline int xlns16_is_positive(xlns16 x) {
    const pbf_t* p = pbf_xlns_ctx();
    uint32_t c = x16_to_pbf(x);
    return (c >= p->mid) && (c != p->max_code);
}

/// Return a sign word: 0x0000 if non-negative, 0x8000 if negative.
static inline xlns16 xlns16_sign(xlns16 x) {
    return (xlns16)(xlns16_is_negative(x) ? 0x8000u : 0x0000u);
}

/// Copy sign of y onto magnitude of x.
static inline xlns16 xlns16_copysign(xlns16 x, xlns16 y) {
    xlns16 ax = xlns16_abs(x);
    return xlns16_is_negative(y) ? xlns16_neg(ax) : ax;
}

// ── Comparisons ──────────────────────────────────────────────────────────────

static inline int xlns16_eq(xlns16 a, xlns16 b) { return a == b; }

static inline int xlns16_gt(xlns16 a, xlns16 b) {
    return pbf_cmp(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)) > 0;
}
static inline int xlns16_lt(xlns16 a, xlns16 b) {
    return pbf_cmp(pbf_xlns_ctx(), x16_to_pbf(a), x16_to_pbf(b)) < 0;
}
static inline int xlns16_ge(xlns16 a, xlns16 b) { return !xlns16_lt(a, b); }
static inline int xlns16_le(xlns16 a, xlns16 b) { return !xlns16_gt(a, b); }

static inline xlns16 xlns16_max(xlns16 a, xlns16 b) { return xlns16_gt(a, b) ? a : b; }
static inline xlns16 xlns16_min(xlns16 a, xlns16 b) { return xlns16_lt(a, b) ? a : b; }

/// Canonical total ordering: larger canon = numerically larger value.
static inline xlns16 xlns16_canon(xlns16 x) { return x; }  // PBF codes are already ordered

// ── Transcendentals (float-boundary implementations) ─────────────────────────

static inline xlns16 xlns16_exp(xlns16 x)  { return fp2xlns16(expf(xlns162fp(x)));  }
static inline xlns16 xlns16_log(xlns16 x)  {
    float f = xlns162fp(x);
    return (f <= 0.0f) ? xlns16_zero : fp2xlns16(logf(f));
}
static inline xlns16 xlns16_exp2(xlns16 x) { return fp2xlns16(exp2f(xlns162fp(x))); }
static inline xlns16 xlns16_log2(xlns16 x) {
    float f = xlns162fp(x);
    return (f <= 0.0f) ? xlns16_zero : fp2xlns16(log2f(f));
}
static inline xlns16 xlns16_pow(xlns16 b, xlns16 e) {
    float fb = xlns162fp(b), fe = xlns162fp(e);
    return (fb <= 0.0f) ? xlns16_zero : fp2xlns16(powf(fb, fe));
}

// ── Activations ──────────────────────────────────────────────────────────────

static inline xlns16 xlns16_relu(xlns16 x) {
    return xlns16_is_negative(x) ? xlns16_zero : x;
}
static inline xlns16 xlns16_sigmoid(xlns16 x) {
    float f = xlns162fp(x);
    return fp2xlns16(1.0f / (1.0f + expf(-f)));
}
static inline xlns16 xlns16_tanh(xlns16 x) { return fp2xlns16(tanhf(xlns162fp(x))); }
static inline xlns16 xlns16_silu(xlns16 x) {
    float f = xlns162fp(x);
    return fp2xlns16(f / (1.0f + expf(-f)));
}
static inline xlns16 xlns16_gelu(xlns16 x) {
    float f = xlns162fp(x);
    float inner = 0.7978845608f * (f + 0.044715f * f * f * f);
    return fp2xlns16(0.5f * f * (1.0f + tanhf(inner)));
}

// ── Vector & array ops ───────────────────────────────────────────────────────

static inline xlns16 xlns16_sum(const xlns16* a, size_t n) {
    if (n == 0) return xlns16_zero;
    xlns16 s = a[0];
    for (size_t i = 1; i < n; i++) s = xlns16_add(s, a[i]);
    return s;
}

static inline xlns16 xlns16_vec_dot(const xlns16* a, const xlns16* b, size_t n) {
    if (n == 0) return xlns16_zero;
    xlns16 s = xlns16_mul(a[0], b[0]);
    for (size_t i = 1; i < n; i++) s = xlns16_add(s, xlns16_mul(a[i], b[i]));
    return s;
}

static inline xlns16 xlns16_max_array(const xlns16* a, size_t n) {
    if (n == 0) return xlns16_zero;
    xlns16 m = a[0];
    for (size_t i = 1; i < n; i++) if (xlns16_gt(a[i], m)) m = a[i];
    return m;
}

static inline xlns16 xlns16_min_array(const xlns16* a, size_t n) {
    if (n == 0) return xlns16_zero;
    xlns16 m = a[0];
    for (size_t i = 1; i < n; i++) if (xlns16_lt(a[i], m)) m = a[i];
    return m;
}

// ── Batch wrappers (forward to pbf_batch_*) ──────────────────────────────────

#define PBF_XLNS_BATCH1(op) \
    for (size_t i = 0; i < n; i++) c[i] = xlns16_##op(a[i])

#define PBF_XLNS_BATCH2(op) \
    for (size_t i = 0; i < n; i++) c[i] = xlns16_##op(a[i], b[i])

static inline void xlns16_batch_add(const xlns16* a, const xlns16* b, xlns16* c, size_t n) { PBF_XLNS_BATCH2(add); }
static inline void xlns16_batch_sub(const xlns16* a, const xlns16* b, xlns16* c, size_t n) { PBF_XLNS_BATCH2(sub); }
static inline void xlns16_batch_mul(const xlns16* a, const xlns16* b, xlns16* c, size_t n) { PBF_XLNS_BATCH2(mul); }
static inline void xlns16_batch_div(const xlns16* a, const xlns16* b, xlns16* c, size_t n) { PBF_XLNS_BATCH2(div); }
static inline void xlns16_batch_neg(const xlns16* a, xlns16* c, size_t n)                  { PBF_XLNS_BATCH1(neg); }
static inline void xlns16_batch_abs(const xlns16* a, xlns16* c, size_t n)                  { PBF_XLNS_BATCH1(abs); }
static inline void xlns16_batch_relu(const xlns16* a, xlns16* c, size_t n)                 { PBF_XLNS_BATCH1(relu); }
static inline void xlns16_batch_sigmoid(const xlns16* a, xlns16* c, size_t n)              { PBF_XLNS_BATCH1(sigmoid); }
static inline void xlns16_batch_tanh(const xlns16* a, xlns16* c, size_t n)                 { PBF_XLNS_BATCH1(tanh); }
static inline void xlns16_batch_silu(const xlns16* a, xlns16* c, size_t n)                 { PBF_XLNS_BATCH1(silu); }
static inline void xlns16_batch_gelu(const xlns16* a, xlns16* c, size_t n)                 { PBF_XLNS_BATCH1(gelu); }

static inline void xlns16_batch_scale(const xlns16* a, xlns16 s, xlns16* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = xlns16_mul(a[i], s);
}

static inline void xlns16_batch_from_float(const float* src, xlns16* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = fp2xlns16(src[i]);
}
static inline void xlns16_batch_to_float(const xlns16* src, float* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = xlns162fp(src[i]);
}

// ── Softmax ──────────────────────────────────────────────────────────────────

static inline void xlns16_softmax_exp(const xlns16* a, xlns16* c, size_t n) {
    xlns16 mx = xlns16_max_array(a, n);
    for (size_t i = 0; i < n; i++) {
        float d = xlns162fp(a[i]) - xlns162fp(mx);
        c[i] = fp2xlns16(expf(d));
    }
}

static inline void xlns16_softmax(const xlns16* a, xlns16* c, size_t n) {
    if (n == 0) return;
    xlns16 mx = xlns16_max_array(a, n);
    for (size_t i = 0; i < n; i++) c[i] = xlns16_exp(xlns16_sub(a[i], mx));
    xlns16 total = xlns16_sum(c, n);
    for (size_t i = 0; i < n; i++) c[i] = xlns16_div(c[i], total);
}

// ── Layer normalization ──────────────────────────────────────────────────────

static inline void xlns16_layernorm(const xlns16* x, xlns16* out,
                                    const xlns16* gamma, const xlns16* beta,
                                    size_t n, float eps) {
    xlns16 mean = xlns16_sum(x, n);
    mean = xlns16_div(mean, fp2xlns16((float)n));
    xlns16 var = xlns16_zero;
    for (size_t i = 0; i < n; i++) {
        xlns16 d = xlns16_sub(x[i], mean);
        var = xlns16_add(var, xlns16_mul(d, d));
    }
    var = xlns16_div(var, fp2xlns16((float)n));
    xlns16 inv_std = fp2xlns16(1.0f / sqrtf(xlns162fp(var) + eps));
    for (size_t i = 0; i < n; i++) {
        out[i] = xlns16_mul(xlns16_sub(x[i], mean), inv_std);
        if (gamma) out[i] = xlns16_mul(out[i], gamma[i]);
        if (beta)  out[i] = xlns16_add(out[i], beta[i]);
    }
}
