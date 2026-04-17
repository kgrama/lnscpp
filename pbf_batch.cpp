// pbf_batch.cpp — Batch, vector, and reduction primitives over PBF arrays.
// Builds on pbf.cpp. All functions take a `const pbf_t*` context.

#pragma once

#include "pbf.cpp"
#include <cstddef>
#include <cmath>

// ── Batch conversions ────────────────────────────────────────────────────────

/// Encode an array of floats to PBF codes.
static inline void pbf_batch_encode(const pbf_t* p, const double* src,
                                    uint32_t* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = pbf_encode(p, src[i]);
}

/// Decode an array of PBF codes to floats.
static inline void pbf_batch_decode(const pbf_t* p, const uint32_t* src,
                                    double* dst, size_t n) {
    for (size_t i = 0; i < n; i++) dst[i] = pbf_decode(p, src[i]);
}

// ── Element-wise arithmetic ──────────────────────────────────────────────────

/// c[i] = a[i] + b[i]
static inline void pbf_batch_add(const pbf_t* p, const uint32_t* a,
                                 const uint32_t* b, uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_add(p, a[i], b[i]);
}

/// c[i] = a[i] - b[i]
static inline void pbf_batch_sub(const pbf_t* p, const uint32_t* a,
                                 const uint32_t* b, uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_sub(p, a[i], b[i]);
}

/// c[i] = a[i] * b[i]
static inline void pbf_batch_mul(const pbf_t* p, const uint32_t* a,
                                 const uint32_t* b, uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_mul(p, a[i], b[i]);
}

/// c[i] = a[i] / b[i]
static inline void pbf_batch_div(const pbf_t* p, const uint32_t* a,
                                 const uint32_t* b, uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_div(p, a[i], b[i]);
}

/// c[i] = a[i] * scalar
static inline void pbf_batch_scale(const pbf_t* p, const uint32_t* a,
                                   uint32_t scalar, uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_mul(p, a[i], scalar);
}

/// c[i] = -a[i]
static inline void pbf_batch_neg(const pbf_t* p, const uint32_t* a,
                                 uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) c[i] = pbf_neg(p, a[i]);
}

/// c[i] = |a[i]|
static inline void pbf_batch_abs(const pbf_t* p, const uint32_t* a,
                                 uint32_t* c, size_t n) {
    for (size_t i = 0; i < n; i++) {
        c[i] = (a[i] != 0 && a[i] < p->mid) ? pbf_neg(p, a[i]) : a[i];
    }
}

// ── Reductions ───────────────────────────────────────────────────────────────

/// Sum of array. Returns 0-code for empty input.
static inline uint32_t pbf_sum(const pbf_t* p, const uint32_t* a, size_t n) {
    if (n == 0) return 0;
    uint32_t s = a[0];
    for (size_t i = 1; i < n; i++) s = pbf_add(p, s, a[i]);
    return s;
}

/// Dot product Σ a[i]*b[i].
static inline uint32_t pbf_dot(const pbf_t* p, const uint32_t* a,
                               const uint32_t* b, size_t n) {
    if (n == 0) return 0;
    uint32_t s = pbf_mul(p, a[0], b[0]);
    for (size_t i = 1; i < n; i++) {
        s = pbf_add(p, s, pbf_mul(p, a[i], b[i]));
    }
    return s;
}

/// Dot product with double inputs (encodes on the fly).
static inline double pbf_vec_dot_f64(const pbf_t* p, const double* a,
                                     const double* b, size_t n) {
    if (n == 0) return 0.0;
    uint32_t s = pbf_mul(p, pbf_encode(p, a[0]), pbf_encode(p, b[0]));
    for (size_t i = 1; i < n; i++) {
        s = pbf_add(p, s, pbf_mul(p, pbf_encode(p, a[i]), pbf_encode(p, b[i])));
    }
    return pbf_decode(p, s);
}

/// Compare two codes: returns -1, 0, +1 in terms of decoded value order.
static inline int pbf_cmp(const pbf_t* p, uint32_t a, uint32_t b) {
    double da = pbf_decode(p, a);
    double db = pbf_decode(p, b);
    return (da < db) ? -1 : (da > db) ? +1 : 0;
}

/// Maximum element. Returns 0-code for empty input.
static inline uint32_t pbf_max(const pbf_t* p, const uint32_t* a, size_t n) {
    if (n == 0) return 0;
    uint32_t m = a[0];
    for (size_t i = 1; i < n; i++) {
        if (pbf_cmp(p, a[i], m) > 0) m = a[i];
    }
    return m;
}

/// Minimum element. Returns 0-code for empty input.
static inline uint32_t pbf_min(const pbf_t* p, const uint32_t* a, size_t n) {
    if (n == 0) return 0;
    uint32_t m = a[0];
    for (size_t i = 1; i < n; i++) {
        if (pbf_cmp(p, a[i], m) < 0) m = a[i];
    }
    return m;
}
