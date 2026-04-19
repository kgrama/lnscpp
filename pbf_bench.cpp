// pbf_bench.cpp — speed benchmark: PBF16 / Tapered PBF16 / xlns32.
//
// Measures ns/op for add, sub, mul, div over pre-encoded operand arrays.
// Operand encoding cost is excluded (paid once, before the timed loop).
// Volatile sinks prevent the loop body from being elided.

#include "pbf_taper.cpp"        // pulls in pbf.cpp
#include "xlns32.cpp"
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <random>
#include <vector>

static const int    N    = 1 << 14;     // operands per array (16384)
static const int    REPS = 200;         // outer reps for stable timing
using clk = std::chrono::steady_clock;

static inline double ns_per_op(clk::duration d, long long ops) {
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(d).count()
         / (double)ops;
}

template <typename Op>
static double bench(Op op) {
    auto t0 = clk::now();
    for (int r = 0; r < REPS; r++) op();
    return ns_per_op(clk::now() - t0, (long long)REPS * N);
}

int main() {
    // Sample operands in the well-defined range of every format.
    std::mt19937 rng(0xBEEF);
    std::uniform_real_distribution<double> dl(-2.0, 2.0);    // 1e-2 .. 1e2
    std::uniform_int_distribution<int>      sg(0, 1);

    std::vector<double>   raw_a(N), raw_b(N);
    for (int i = 0; i < N; i++) {
        raw_a[i] = std::pow(10.0, dl(rng)) * (sg(rng) ? 1.0 : -1.0);
        raw_b[i] = std::pow(10.0, dl(rng)) * (sg(rng) ? 1.0 : -1.0);
    }

    pbf_t       pbf16; pbf16_init(&pbf16);
    pbf_taper_t tap16; pbf_taper16_init(&tap16);

    std::vector<uint32_t> p_a(N), p_b(N), p_c(N);
    std::vector<uint32_t> t_a(N), t_b(N), t_c(N);
    std::vector<xlns32>   x_a(N), x_b(N), x_c(N);
    for (int i = 0; i < N; i++) {
        p_a[i] = pbf_encode      (&pbf16, raw_a[i]);
        p_b[i] = pbf_encode      (&pbf16, raw_b[i]);
        t_a[i] = pbf_taper_encode(&tap16, raw_a[i]);
        t_b[i] = pbf_taper_encode(&tap16, raw_b[i]);
        x_a[i] = fp2xlns32       ((float)raw_a[i]);
        x_b[i] = fp2xlns32       ((float)raw_b[i]);
    }

    volatile uint32_t sink_u = 0;
    volatile xlns32   sink_x = 0;

    typedef uint32_t (*pbf_fn_t)(const pbf_t*, uint32_t, uint32_t);
    typedef uint32_t (*tap_fn_t)(const pbf_taper_t*, uint32_t, uint32_t);
    typedef xlns32   (*x32_fn_t)(xlns32, xlns32);

    auto t0 = clk::now();
    (void)t0;

    struct { pbf_fn_t fn; const char* name; double ns; } pbf_ops[] = {
        {pbf_add,"add",0}, {pbf_sub,"sub",0}, {pbf_mul,"mul",0}, {pbf_div,"div",0},
    };
    struct { tap_fn_t fn; const char* name; double ns; } tap_ops[] = {
        {pbf_taper_add,"add",0}, {pbf_taper_sub,"sub",0},
        {pbf_taper_mul,"mul",0}, {pbf_taper_div,"div",0},
    };
    struct { x32_fn_t fn; const char* name; double ns; } x32_ops[] = {
        {xlns32_add,"add",0},
        {[](xlns32 a, xlns32 b){ return xlns32_sub(a,b); },"sub",0},
        {xlns32_mul,"mul",0}, {xlns32_div,"div",0},
    };

    for (auto& op : pbf_ops) {
        pbf_fn_t fn = op.fn;
        op.ns = bench([&]{
            for (int i = 0; i < N; i++) p_c[i] = fn(&pbf16, p_a[i], p_b[i]);
            sink_u = p_c[N-1];
        });
    }
    for (auto& op : tap_ops) {
        tap_fn_t fn = op.fn;
        op.ns = bench([&]{
            for (int i = 0; i < N; i++) t_c[i] = fn(&tap16, t_a[i], t_b[i]);
            sink_u = t_c[N-1];
        });
    }
    for (auto& op : x32_ops) {
        x32_fn_t fn = op.fn;
        op.ns = bench([&]{
            for (int i = 0; i < N; i++) x_c[i] = fn(x_a[i], x_b[i]);
            sink_x = x_c[N-1];
        });
    }

    (void)sink_u; (void)sink_x;
    double pbf_add_ns = pbf_ops[0].ns, pbf_sub_ns = pbf_ops[1].ns;
    double pbf_mul_ns = pbf_ops[2].ns, pbf_div_ns = pbf_ops[3].ns;
    double tap_add_ns = tap_ops[0].ns, tap_sub_ns = tap_ops[1].ns;
    double tap_mul_ns = tap_ops[2].ns, tap_div_ns = tap_ops[3].ns;
    double x32_add_ns = x32_ops[0].ns, x32_sub_ns = x32_ops[1].ns;
    double x32_mul_ns = x32_ops[2].ns, x32_div_ns = x32_ops[3].ns;

    std::printf("=================================================================\n");
    std::printf("  Speed benchmark (ns/op, lower is better)\n");
    std::printf("  N=%d operands, REPS=%d, range 1e-2..1e2\n", N, REPS);
    std::printf("=================================================================\n");
    std::printf("  %-14s   %8s  %8s  %8s  %8s\n", "format", "add", "sub", "mul", "div");
    std::printf("  %-14s   %7.1fn %7.1fn %7.1fn %7.1fn\n",
                "PBF16",         pbf_add_ns, pbf_sub_ns, pbf_mul_ns, pbf_div_ns);
    std::printf("  %-14s   %7.1fn %7.1fn %7.1fn %7.1fn\n",
                "Tapered PBF16", tap_add_ns, tap_sub_ns, tap_mul_ns, tap_div_ns);
    std::printf("  %-14s   %7.1fn %7.1fn %7.1fn %7.1fn\n",
                "xlns32 (LUT)",  x32_add_ns, x32_sub_ns, x32_mul_ns, x32_div_ns);

    std::printf("\n  ratios (PBF16 / xlns32, <1 = PBF faster)\n");
    std::printf("  %-14s   %8.2f  %8.2f  %8.2f  %8.2f\n", "PBF16/xlns32",
                pbf_add_ns/x32_add_ns, pbf_sub_ns/x32_sub_ns,
                pbf_mul_ns/x32_mul_ns, pbf_div_ns/x32_div_ns);
    return 0;
}
