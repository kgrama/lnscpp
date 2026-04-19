// pbf_taper_test.cpp — SNR comparison: flat PBF16 vs tapered PBF16.
//
// Sweeps operands across three regimes:
//   1. Near-unity (|x| ∈ [0.5, 2]) — where tapered should dominate
//   2. Mid-range (|x| ∈ [1e-2, 1e2])
//   3. Wide (|x| ∈ [1e-7, 1e7]) — where tapered should slightly underperform
// and reports SNR for roundtrip / mul / div / add for each format and regime.

#include "pbf_taper.cpp"   // transitively includes pbf.cpp
#include <cstdio>
#include <cmath>
#include <random>

struct snr_acc {
    double sum_sq_signal = 0;
    double sum_sq_error  = 0;
    void add(double ref, double got) {
        if (!std::isfinite(ref) || !std::isfinite(got)) return;
        sum_sq_signal += ref * ref;
        sum_sq_error  += (got - ref) * (got - ref);
    }
    double db() const {
        if (sum_sq_error == 0) return 200.0;
        return 10.0 * std::log10(sum_sq_signal / sum_sq_error);
    }
};

struct regime { const char* name; double lo_log10; double hi_log10; };

static const regime REGIMES[] = {
    { "near-unity [0.5, 2]",     -0.301,  0.301 },
    { "mid-range  [1e-2, 1e2]",  -2.0,    2.0   },
    { "wide       [1e-7, 1e7]",  -7.0,    7.0   },
};

static double rand_signed(std::mt19937& rng, double lo_log10, double hi_log10) {
    std::uniform_real_distribution<double> dl(lo_log10, hi_log10);
    std::uniform_int_distribution<int>      sg(0, 1);
    double v = std::pow(10.0, dl(rng));
    return sg(rng) ? v : -v;
}

template <typename Pbf, typename Encode, typename Decode,
          typename Mul, typename Div, typename Add, typename Sub>
static void run_format(const char* label, const Pbf& p, int n_samples,
                       Encode enc, Decode dec, Mul mul, Div div, Add add, Sub sub) {
    std::printf("\n%s\n", label);
    std::printf("  %-26s   roundtrip      mul      div      add\n", "regime");
    for (const auto& r : REGIMES) {
        std::mt19937 rng(0xC0FFEE);
        snr_acc rt, mu, dv, ad;
        for (int i = 0; i < n_samples; i++) {
            double x  = rand_signed(rng, r.lo_log10, r.hi_log10);
            double y  = rand_signed(rng, r.lo_log10, r.hi_log10);
            // Roundtrip
            rt.add(x, dec(&p, enc(&p, x)));
            // Mul (skip if product overflows the format's outer bound)
            double xy = x * y;
            if (std::fabs(xy) < 1e8 && std::fabs(xy) > 1e-8) {
                mu.add(xy, dec(&p, mul(&p, enc(&p, x), enc(&p, y))));
            }
            // Div
            if (y != 0.0) {
                double xd = x / y;
                if (std::fabs(xd) < 1e8 && std::fabs(xd) > 1e-8) {
                    dv.add(xd, dec(&p, div(&p, enc(&p, x), enc(&p, y))));
                }
            }
            // Add
            ad.add(x + y, dec(&p, add(&p, enc(&p, x), enc(&p, y))));
            (void)sub;
        }
        std::printf("  %-26s   %5.1f dB  %5.1f dB  %5.1f dB  %5.1f dB\n",
                    r.name, rt.db(), mu.db(), dv.db(), ad.db());
    }
}

int main() {
    pbf_t flat8, flat16;
    pbf8_init (&flat8);
    pbf16_init(&flat16);
    pbf_taper_t tap8, tap16;
    pbf_taper8_init (&tap8);
    pbf_taper16_init(&tap16);

    std::printf("=================================================================\n");
    std::printf("  Flat PBF vs Tapered PBF — step ratios T0:T1:T2:T3 = 8:1:1:8\n");
    std::printf("  (denser around |x| = 1.0; W16 boundaries match P32_Taper.v)\n");
    std::printf("=================================================================\n");

    run_format("Flat PBF8    (uniform log spacing, range 1e-4..1e4)", flat8, 4000,
               pbf_encode, pbf_decode, pbf_mul, pbf_div, pbf_add, pbf_sub);
    run_format("Tapered PBF8 (4 tapers, range 1e-4..1e4)",            tap8,  4000,
               pbf_taper_encode, pbf_taper_decode, pbf_taper_mul,
               pbf_taper_div, pbf_taper_add, pbf_taper_sub);

    run_format("Flat PBF16    (uniform log spacing, range 1e-8..1e8)", flat16, 4000,
               pbf_encode, pbf_decode, pbf_mul, pbf_div, pbf_add, pbf_sub);
    run_format("Tapered PBF16 (4 tapers, range 1e-8..1e8)",            tap16,  4000,
               pbf_taper_encode, pbf_taper_decode, pbf_taper_mul,
               pbf_taper_div, pbf_taper_add, pbf_taper_sub);

    return 0;
}
