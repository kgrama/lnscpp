// pbf_wasm.cpp — emscripten entry points for the fractal demo.
//
// Exposes the existing PBF and tapered-PBF arithmetic from pbf.cpp /
// pbf_taper.cpp directly to JavaScript. No reimplementation, no Q13
// compromise — the WASM module runs the same Q32 CF as the C++ tests.
//
// Build:
//   emcc -O3 -std=c++11 pbf_wasm.cpp -o pbf_wasm.js \
//        -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME=createPbfModule \
//        -s EXPORTED_FUNCTIONS='["_pbf_init_w","_pbf_taper_init_w",
//                                "_pbf_encode_w","_pbf_decode_w",
//                                "_pbf_add_w","_pbf_sub_w","_pbf_mul_w","_pbf_div_w",
//                                "_pbf_t_encode_w","_pbf_t_decode_w",
//                                "_pbf_t_add_w","_pbf_t_sub_w","_pbf_t_mul_w","_pbf_t_div_w",
//                                "_render_mandelbrot","_render_julia","_render_burning",
//                                "_malloc","_free"]' \
//        -s EXPORTED_RUNTIME_METHODS='["HEAPU32","HEAPU8","HEAPF64"]'

#include "../pbf_taper.cpp"      // transitively pulls pbf.cpp
#include <cstdint>

extern "C" {

// Persistent format descriptors. JS calls *_init_w once with parameters,
// subsequent calls reuse the descriptor.
static pbf_t       g_flat;
static pbf_taper_t g_taper;

// ── Format setup ────────────────────────────────────────────────────────────
void pbf_init_w(int n_bits, double v_min, double v_max) {
    pbf_init(&g_flat, n_bits, v_min, v_max);
}

void pbf_taper_init_w(int n_bits, double v_max) {
    pbf_taper_init(&g_taper, n_bits, v_max);
}

// ── Flat PBF scalar ops (thin wrappers; JS sees uint32_t codes) ─────────────
uint32_t pbf_encode_w(double v)               { return pbf_encode(&g_flat, v); }
double   pbf_decode_w(uint32_t c)             { return pbf_decode(&g_flat, c); }
uint32_t pbf_add_w   (uint32_t a, uint32_t b) { return pbf_add(&g_flat, a, b); }
uint32_t pbf_sub_w   (uint32_t a, uint32_t b) { return pbf_sub(&g_flat, a, b); }
uint32_t pbf_mul_w   (uint32_t a, uint32_t b) { return pbf_mul(&g_flat, a, b); }
uint32_t pbf_div_w   (uint32_t a, uint32_t b) { return pbf_div(&g_flat, a, b); }

// ── Tapered PBF scalar ops ──────────────────────────────────────────────────
uint32_t pbf_t_encode_w(double v)               { return pbf_taper_encode(&g_taper, v); }
double   pbf_t_decode_w(uint32_t c)             { return pbf_taper_decode(&g_taper, c); }
uint32_t pbf_t_add_w   (uint32_t a, uint32_t b) { return pbf_taper_add(&g_taper, a, b); }
uint32_t pbf_t_sub_w   (uint32_t a, uint32_t b) { return pbf_taper_sub(&g_taper, a, b); }
uint32_t pbf_t_mul_w   (uint32_t a, uint32_t b) { return pbf_taper_mul(&g_taper, a, b); }
uint32_t pbf_t_div_w   (uint32_t a, uint32_t b) { return pbf_taper_div(&g_taper, a, b); }

// ── Whole-image fractal kernels ─────────────────────────────────────────────
// Each kernel writes one byte per pixel (iteration count clamped to 0..255)
// into `out`, which JS pre-allocates and reads back as a Uint8ClampedArray.
// `is_taper` selects the format (0 = flat, 1 = tapered).
//
// The full fractal iteration runs in the WASM module — no per-pixel JS↔WASM
// crossings inside the inner loop, where the per-call overhead would
// dominate over the PBF op itself.

static inline uint8_t iter_byte(int n, int max_iter) {
    if (n >= max_iter) return 0;
    return (uint8_t)((n * 255) / max_iter);
}

// Mandelbrot: z := z² + c, c = pixel
void render_mandelbrot(uint8_t* out, int w, int h,
                       double cx, double cy, double scale,
                       int max_iter, int is_taper) {
    const double aspect = (double)w / (double)h;
    for (int py = 0; py < h; py++) {
        for (int px = 0; px < w; px++) {
            double u = ((double)px / (double)w  - 0.5) * aspect;
            double v = ((double)py / (double)h  - 0.5);
            double cre = cx + u * scale;
            double cim = cy + v * scale;

            uint32_t c_re, c_im, z_re, z_im, four;
            if (is_taper) {
                c_re = pbf_taper_encode(&g_taper, cre);
                c_im = pbf_taper_encode(&g_taper, cim);
                four = pbf_taper_encode(&g_taper, 4.0);
            } else {
                c_re = pbf_encode(&g_flat, cre);
                c_im = pbf_encode(&g_flat, cim);
                four = pbf_encode(&g_flat, 4.0);
            }
            z_re = c_re; z_im = c_im;

            int n = 0;
            for (; n < max_iter; n++) {
                uint32_t r2, i2, ri, mag2;
                if (is_taper) {
                    r2   = pbf_taper_mul(&g_taper, z_re, z_re);
                    i2   = pbf_taper_mul(&g_taper, z_im, z_im);
                    ri   = pbf_taper_mul(&g_taper, z_re, z_im);
                    mag2 = pbf_taper_add(&g_taper, r2, i2);
                    if (mag2 > four) break;     // SBP code compare = magnitude compare for positives
                    z_re = pbf_taper_add(&g_taper, pbf_taper_sub(&g_taper, r2, i2), c_re);
                    z_im = pbf_taper_add(&g_taper, pbf_taper_add(&g_taper, ri, ri), c_im);
                } else {
                    r2   = pbf_mul(&g_flat, z_re, z_re);
                    i2   = pbf_mul(&g_flat, z_im, z_im);
                    ri   = pbf_mul(&g_flat, z_re, z_im);
                    mag2 = pbf_add(&g_flat, r2, i2);
                    if (mag2 > four) break;
                    z_re = pbf_add(&g_flat, pbf_sub(&g_flat, r2, i2), c_re);
                    z_im = pbf_add(&g_flat, pbf_add(&g_flat, ri, ri), c_im);
                }
            }
            out[py * w + px] = iter_byte(n, max_iter);
        }
    }
}

// Julia: z := z² + c, c = fixed seed
void render_julia(uint8_t* out, int w, int h,
                  double cx, double cy, double scale,
                  double seed_re, double seed_im,
                  int max_iter, int is_taper) {
    const double aspect = (double)w / (double)h;
    uint32_t c_re, c_im, four;
    if (is_taper) {
        c_re = pbf_taper_encode(&g_taper, seed_re);
        c_im = pbf_taper_encode(&g_taper, seed_im);
        four = pbf_taper_encode(&g_taper, 4.0);
    } else {
        c_re = pbf_encode(&g_flat, seed_re);
        c_im = pbf_encode(&g_flat, seed_im);
        four = pbf_encode(&g_flat, 4.0);
    }
    for (int py = 0; py < h; py++) {
        for (int px = 0; px < w; px++) {
            double u = ((double)px / (double)w  - 0.5) * aspect;
            double v = ((double)py / (double)h  - 0.5);
            double zre = cx + u * scale;
            double zim = cy + v * scale;
            uint32_t z_re = is_taper ? pbf_taper_encode(&g_taper, zre)
                                     : pbf_encode(&g_flat, zre);
            uint32_t z_im = is_taper ? pbf_taper_encode(&g_taper, zim)
                                     : pbf_encode(&g_flat, zim);
            int n = 0;
            for (; n < max_iter; n++) {
                uint32_t r2, i2, ri, mag2;
                if (is_taper) {
                    r2   = pbf_taper_mul(&g_taper, z_re, z_re);
                    i2   = pbf_taper_mul(&g_taper, z_im, z_im);
                    ri   = pbf_taper_mul(&g_taper, z_re, z_im);
                    mag2 = pbf_taper_add(&g_taper, r2, i2);
                    if (mag2 > four) break;
                    z_re = pbf_taper_add(&g_taper, pbf_taper_sub(&g_taper, r2, i2), c_re);
                    z_im = pbf_taper_add(&g_taper, pbf_taper_add(&g_taper, ri, ri), c_im);
                } else {
                    r2   = pbf_mul(&g_flat, z_re, z_re);
                    i2   = pbf_mul(&g_flat, z_im, z_im);
                    ri   = pbf_mul(&g_flat, z_re, z_im);
                    mag2 = pbf_add(&g_flat, r2, i2);
                    if (mag2 > four) break;
                    z_re = pbf_add(&g_flat, pbf_sub(&g_flat, r2, i2), c_re);
                    z_im = pbf_add(&g_flat, pbf_add(&g_flat, ri, ri), c_im);
                }
            }
            out[py * w + px] = iter_byte(n, max_iter);
        }
    }
}

// Burning Ship: z := (|re| + i|im|)² + c
void render_burning(uint8_t* out, int w, int h,
                    double cx, double cy, double scale,
                    int max_iter, int is_taper) {
    const double aspect = (double)w / (double)h;
    for (int py = 0; py < h; py++) {
        for (int px = 0; px < w; px++) {
            double u = ((double)px / (double)w  - 0.5) * aspect;
            double v = ((double)py / (double)h  - 0.5);
            double cre = cx + u * scale;
            double cim = cy + v * scale;
            uint32_t c_re, c_im, z_re, z_im, four;
            if (is_taper) {
                c_re = pbf_taper_encode(&g_taper, cre);
                c_im = pbf_taper_encode(&g_taper, cim);
                four = pbf_taper_encode(&g_taper, 4.0);
            } else {
                c_re = pbf_encode(&g_flat, cre);
                c_im = pbf_encode(&g_flat, cim);
                four = pbf_encode(&g_flat, 4.0);
            }
            z_re = c_re; z_im = c_im;

            int n = 0;
            for (; n < max_iter; n++) {
                // |z| = (|re|, |im|) — strip sign by setting code into the
                // positive half (codes >= mid encode positives in SBP).
                uint32_t azr = is_taper
                    ? (z_re < g_taper.mid ? (g_taper.max_code - z_re) : z_re)
                    : (z_re < g_flat.mid  ? (g_flat.max_code  - z_re) : z_re);
                uint32_t azi = is_taper
                    ? (z_im < g_taper.mid ? (g_taper.max_code - z_im) : z_im)
                    : (z_im < g_flat.mid  ? (g_flat.max_code  - z_im) : z_im);

                uint32_t r2, i2, ri, mag2;
                if (is_taper) {
                    r2   = pbf_taper_mul(&g_taper, azr, azr);
                    i2   = pbf_taper_mul(&g_taper, azi, azi);
                    ri   = pbf_taper_mul(&g_taper, azr, azi);
                    mag2 = pbf_taper_add(&g_taper, r2, i2);
                    if (mag2 > four) break;
                    z_re = pbf_taper_add(&g_taper, pbf_taper_sub(&g_taper, r2, i2), c_re);
                    z_im = pbf_taper_add(&g_taper, pbf_taper_add(&g_taper, ri, ri), c_im);
                } else {
                    r2   = pbf_mul(&g_flat, azr, azr);
                    i2   = pbf_mul(&g_flat, azi, azi);
                    ri   = pbf_mul(&g_flat, azr, azi);
                    mag2 = pbf_add(&g_flat, r2, i2);
                    if (mag2 > four) break;
                    z_re = pbf_add(&g_flat, pbf_sub(&g_flat, r2, i2), c_re);
                    z_im = pbf_add(&g_flat, pbf_add(&g_flat, ri, ri), c_im);
                }
            }
            out[py * w + px] = iter_byte(n, max_iter);
        }
    }
}

}  // extern "C"
