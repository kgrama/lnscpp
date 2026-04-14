#ifdef USE_SBP
#  include "../pbf16_sbp.cpp"
#  define IMPL "SBP"
#elif defined(USE_BFLOAT)
#  define IMPL "BFLOAT16"
#else
#  include "../pbf16.cpp"
#  define IMPL "PBF"
#endif

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <ctime>

// ── bfloat16: truncate float32 to top 16 bits ────────────────────────────────
static uint16_t f2bf(float f) {
    uint32_t u; __builtin_memcpy(&u, &f, 4);
    return (uint16_t)(u >> 16);
}
static float bf2f(uint16_t b) {
    uint32_t u = (uint32_t)b << 16;
    float f; __builtin_memcpy(&f, &u, 4);
    return f;
}
static uint16_t bf_add(uint16_t a, uint16_t b) { return f2bf(bf2f(a) + bf2f(b)); }

// ── dispatch: same call site for all three impls ──────────────────────────────
static uint16_t enc(float x) {
#ifdef USE_BFLOAT
    return f2bf(x);
#else
    return fp2xlns16(x);
#endif
}
static float dec(uint16_t x) {
#ifdef USE_BFLOAT
    return bf2f(x);
#else
    return xlns162fp(x);
#endif
}
static uint16_t do_add(uint16_t a, uint16_t b) {
#ifdef USE_BFLOAT
    return bf_add(a, b);
#else
    return xlns16_add(a, b);
#endif
}

// ── RNG ───────────────────────────────────────────────────────────────────────
static uint32_t lcg_s = 42;
static float rnd(float lo, float hi) {
    lcg_s = lcg_s * 1664525u + 1013904223u;
    return lo + ((float)(lcg_s >> 8) / (float)(1<<24)) * (hi - lo);
}

static int pass = 0, fail = 0;
static void chk(const char* lbl, float got, float ref, float tol = 2.5f) {
    bool nz = fabsf(ref) < 1e-5f;
    float e  = nz ? fabsf(got)*100.f : fabsf(got-ref)/fabsf(ref)*100.f;
    bool ok  = nz ? fabsf(got) < 5e-3f : e <= tol;
    ok ? pass++ : fail++;
    if (nz) printf("  [%s] %-28s  ref=%9.4f  got=%9.4f  err=~0\n",   ok?"OK  ":"FAIL",lbl,ref,got);
    else    printf("  [%s] %-28s  ref=%9.4f  got=%9.4f  err=%.2f%%\n",ok?"OK  ":"FAIL",lbl,ref,got,e);
}

static void snr_op(const char* nm, int mode) {
    double sig=0, noise=0; int n=0;
    lcg_s = 0xdeadbeef;
    for (int i = 0; i < 20000; i++) {
        float a = rnd(-16.f, 16.f), b = rnd(-16.f, 16.f);
        if (fabsf(a)<1e-3f || fabsf(b)<1e-3f) continue;
        if (mode==3 && fabsf(b)<0.05f) continue;
        float ref; uint16_t xa=enc(a), xb=enc(b), xr;
#ifndef USE_BFLOAT
        switch(mode){
            case 0: ref=a+b; xr=xlns16_add(xa,xb); break;
            case 1: ref=a-b; xr=xlns16_sub(xa,xb); break;
            case 2: ref=a*b; xr=xlns16_mul(xa,xb); break;
            default:ref=a/b; xr=xlns16_div(xa,xb); break;
        }
#else
        switch(mode){
            case 0: ref=a+b; xr=f2bf(bf2f(xa)+bf2f(xb)); break;
            case 1: ref=a-b; xr=f2bf(bf2f(xa)-bf2f(xb)); break;
            case 2: ref=a*b; xr=f2bf(bf2f(xa)*bf2f(xb)); break;
            default:ref=a/b; xr=f2bf(bf2f(xa)/bf2f(xb)); break;
        }
#endif
        if (!std::isfinite(ref) || fabsf(ref)>500.f) continue;
        float got = dec(xr);
        sig   += (double)ref*ref + 1e-8;
        noise += (double)(got-ref)*(got-ref);
        n++;
    }
    printf("  %-5s  SNR: %6.1f dB  (n=%d)\n", nm, 10*log10(sig/(noise+1e-300)), n);
}

static void test_arith() {
    printf("\n--- Basic arithmetic ---\n");
    float cs[][2]={{1,1},{3,4},{-3,4},{0.5f,0.5f},{7,0.125f},{-5,3},{8,0.5f}};
    for (auto& c : cs) {
        char l[24];
        snprintf(l,24,"%g+%g",c[0],c[1]);
        chk(l, dec(do_add(enc(c[0]),enc(c[1]))), c[0]+c[1],
            fabsf(c[0]+c[1])<0.1f?50.f:2.5f);
    }
}

static void test_snr() {
    printf("\n--- SNR vs IEEE float (N=20000, range ±16) ---\n");
    snr_op("add",0); snr_op("sub",1); snr_op("mul",2); snr_op("div",3);
}

static void test_roundtrip() {
    printf("\n--- Round-trip SNR ---\n");
    float vs[]={0.001f,.01f,.1f,.5f,1,2,4,8,32,-1,-2,-8,1.5f,3.14f,0.707f};
    double ps=0,rs=0; int N=sizeof(vs)/sizeof(vs[0]);
    for(float v:vs){float g=dec(enc(v));ps+=(double)(g-v)*(g-v);rs+=(double)v*v;}
    printf("  SNR: %.1f dB\n", 10*log10(rs/N/(ps/N+1e-300)));
}

// ── cancellation: the main comparison ────────────────────────────────────────
static void test_cancellation() {
    printf("\n--- Near-cancellation: a + (-a*(1+eps)), eps in [0.001, 0.5] ---\n");
    printf("  %-10s  %-10s  %-10s  %-10s  %-8s\n","a","b","ref","got","abserr");

    double sig=0, noise=0; int n=0;
    lcg_s = 0xabcd1234;

    // print a few samples
    int samples = 0;
    uint32_t saved = lcg_s;
    for (int i = 0; i < 20; i++) {
        float a   = rnd(0.5f, 8.f);
        float eps = rnd(0.001f, 0.5f);
        float b   = -a * (1.f + eps);
        float ref = a + b;
        float got = dec(do_add(enc(a), enc(b)));
        printf("  %-10.4f  %-10.4f  %-10.5f  %-10.5f  %.5f\n",
               a, b, ref, got, fabsf(got-ref));
    }

    // full SNR sweep
    lcg_s = saved;
    for (int i = 0; i < 10000; i++) {
        float a   = rnd(0.5f, 8.f);
        float eps = rnd(0.001f, 0.5f);
        float b   = -a * (1.f + eps);
        float ref = a + b;
        float got = dec(do_add(enc(a), enc(b)));
        sig   += (double)ref*ref + 1e-8;
        noise += (double)(got-ref)*(got-ref);
        n++;
    }
    printf("  Cancellation SNR: %.1f dB  (n=%d)\n",
           10*log10(sig/(noise+1e-300)), n);
}

int main(){
    printf("==========================================================\n");
    printf("  %s\n", IMPL);
    printf("==========================================================\n");
    test_arith();
    test_snr();
    test_roundtrip();
    test_cancellation();
    printf("\n  SUMMARY: %d passed, %d failed\n", pass, fail);
    return fail ? 1 : 0;
}
