// pbftest.cpp — SNR and CF-accuracy demo for pbf.cpp.
// Mirrors the Python __main__ block of pbf.py.
//
// Build: g++ -std=c++11 -O2 pbftest.cpp -o pbftest && ./pbftest

#include "pbf.cpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#include <algorithm>

// ── Deterministic PRNG helpers (std::mt19937) ────────────────────────────────

struct pbf_rng {
    std::mt19937_64 eng;
    pbf_rng(uint64_t seed) : eng(seed) {}
    double uniform(double lo, double hi) {
        std::uniform_real_distribution<double> d(lo, hi);
        return d(eng);
    }
    int coin() {
        std::uniform_int_distribution<int> d(0, 1);
        return d(eng) ? +1 : -1;
    }
};

// ── SNR accumulator ──────────────────────────────────────────────────────────

struct snr_acc {
    double signal_sq = 0.0;
    double noise_sq  = 0.0;
    void add(double expected, double got) {
        signal_sq += expected * expected;
        double e = got - expected;
        noise_sq += e * e;
    }
    double db() const {
        double n = noise_sq > 1e-100 ? noise_sq : 1e-100;
        return 10.0 * std::log10(signal_sq / n);
    }
};

// ── SNR suite: roundtrip, add, sub, mul, div ─────────────────────────────────

struct snr_results { double roundtrip, add, sub, mul, div; };

static snr_results snr_test(const pbf_t& p, int n_samples, uint64_t seed) {
    pbf_rng rng(seed);
    snr_results r{};

    double log_min = std::log10(p.v_min);
    double log_max = std::log10(p.v_max);
    double log_range = log_max - log_min;
    double safe_min = log_min + log_range * 0.1;
    double safe_max = log_max - log_range * 0.1;
    double mid       = (safe_min + safe_max) / 2.0;
    double add_range = (safe_max - safe_min) / 3.0;
    double mul_range = std::min(add_range / 2.0, log_range / 4.0);

    // Roundtrip
    snr_acc a;
    for (int i = 0; i < n_samples; i++) {
        double v = std::pow(10.0, rng.uniform(safe_min, safe_max)) * rng.coin();
        a.add(v, pbf_decode(&p, pbf_encode(&p, v)));
    }
    r.roundtrip = a.db();

    // Addition
    a = snr_acc();
    for (int i = 0; i < n_samples; i++) {
        double x = std::pow(10.0, rng.uniform(mid - add_range, mid + add_range));
        double y = std::pow(10.0, rng.uniform(mid - add_range, mid + add_range));
        double expected = x + y;
        double got = pbf_decode(&p, pbf_add(&p, pbf_encode(&p, x), pbf_encode(&p, y)));
        a.add(expected, got);
    }
    r.add = a.db();

    // Subtraction: a > b, both positive
    a = snr_acc();
    for (int i = 0; i < n_samples; i++) {
        double x = std::pow(10.0, rng.uniform(mid, mid + add_range));
        double y = std::pow(10.0, rng.uniform(mid - add_range, mid));
        if (x <= y) { double t = y * 1.1; y = x; x = t; }
        double expected = x - y;
        double got = pbf_decode(&p, pbf_sub(&p, pbf_encode(&p, x), pbf_encode(&p, y)));
        a.add(expected, got);
    }
    r.sub = a.db();

    // Multiplication
    a = snr_acc();
    for (int i = 0; i < n_samples; i++) {
        double x = std::pow(10.0, rng.uniform(mid - mul_range, mid + mul_range));
        double y = std::pow(10.0, rng.uniform(mid - mul_range, mid + mul_range));
        double expected = x * y;
        double got = pbf_decode(&p, pbf_mul(&p, pbf_encode(&p, x), pbf_encode(&p, y)));
        a.add(expected, got);
    }
    r.mul = a.db();

    // Division
    a = snr_acc();
    for (int i = 0; i < n_samples; i++) {
        double x = std::pow(10.0, rng.uniform(mid - mul_range, mid + mul_range));
        double y = std::pow(10.0, rng.uniform(mid - mul_range, mid + mul_range));
        double expected = x / y;
        double got = pbf_decode(&p, pbf_div(&p, pbf_encode(&p, x), pbf_encode(&p, y)));
        a.add(expected, got);
    }
    r.div = a.db();

    return r;
}

// ── Error distribution for addition ──────────────────────────────────────────

struct err_stats {
    double min, p5, p25, median, p75, p95, max, mean, rms;
};

static err_stats error_distribution(const pbf_t& p, int n_samples, uint64_t seed) {
    pbf_rng rng(seed);
    double log_min = std::log10(p.v_min) + 1.0;
    double log_max = std::log10(p.v_max) - 1.0;

    std::vector<double> errors;
    errors.reserve(n_samples);
    for (int i = 0; i < n_samples; i++) {
        double x = std::pow(10.0, rng.uniform(log_min, log_max));
        double y = std::pow(10.0, rng.uniform(log_min, log_max));
        double expected = x + y;
        double got = pbf_decode(&p, pbf_add(&p, pbf_encode(&p, x), pbf_encode(&p, y)));
        if (expected != 0.0) errors.push_back((got - expected) / expected * 100.0);
    }
    std::sort(errors.begin(), errors.end());

    err_stats s{};
    int n = (int)errors.size();
    auto pct = [&](double q) { return errors[std::min(n - 1, (int)(n * q))]; };
    s.min    = errors.front();
    s.p5     = pct(0.05);
    s.p25    = pct(0.25);
    s.median = errors[n / 2];
    s.p75    = pct(0.75);
    s.p95    = pct(0.95);
    s.max    = errors.back();
    double sum = 0.0, sumsq = 0.0;
    for (double e : errors) { sum += e; sumsq += e * e; }
    s.mean = sum / n;
    s.rms  = std::sqrt(sumsq / n);
    return s;
}

// ── Main: replicate the pbf.py demo ──────────────────────────────────────────

int main() {
    std::printf("=================================================================\n");
    std::printf("  PBF: Table-Free Continued Fraction Edition  (C++ port)\n");
    std::printf("=================================================================\n");

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("SNR by Format (higher = better, ~6 dB per bit)\n");
    std::printf("-----------------------------------------------------------------\n");

    pbf_t p8, p16;
    pbf8_init(&p8);
    pbf16_init(&p16);

    struct entry { const char* name; const pbf_t* p; } fmts[] = {
        { "PBF8 ", &p8 },
        { "PBF16", &p16 },
    };

    std::printf("\n%-8s %12s %12s %12s %12s %12s\n",
                "Format", "Roundtrip", "Add", "Sub", "Mul", "Div");
    std::printf("-----------------------------------------------------------------\n");

    for (auto& e : fmts) {
        snr_results r = snr_test(*e.p, 2000, 42);
        std::printf("%-8s", e.name);
        const double ops[] = { r.roundtrip, r.add, r.sub, r.mul, r.div };
        for (double db : ops) {
            double bits = db / 6.02;
            std::printf(" %5.1fdB/%.1fb", db, bits);
        }
        std::printf("\n");
    }

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("Error Distribution for Addition (PBF8)\n");
    std::printf("-----------------------------------------------------------------\n\n");

    err_stats d = error_distribution(p8, 5000, 42);
    std::printf("  Min:    %+7.2f%%\n", d.min);
    std::printf("  5th %%:  %+7.2f%%\n", d.p5);
    std::printf("  25th %%: %+7.2f%%\n", d.p25);
    std::printf("  Median: %+7.2f%%\n", d.median);
    std::printf("  75th %%: %+7.2f%%\n", d.p75);
    std::printf("  95th %%: %+7.2f%%\n", d.p95);
    std::printf("  Max:    %+7.2f%%\n", d.max);
    std::printf("  Mean:   %+7.2f%%\n", d.mean);
    std::printf("  RMS:    %7.2f%%\n",  d.rms);

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf("CF Primitive Accuracy\n");
    std::printf("-----------------------------------------------------------------\n");

    struct test { const char* name; double got; double exact; };
    double q = (double)PBF_Q;
    test tests[] = {
        { "ln(2)",  pbf_ln_cf(2 * PBF_Q)                           / q, std::log(2.0)  },
        { "ln(e)",  pbf_ln_cf((int64_t)(std::exp(1.0) * q))        / q, 1.0            },
        { "ln(10)", pbf_ln_cf(10 * PBF_Q)                          / q, std::log(10.0) },
        { "exp(1)", pbf_exp_cf(PBF_Q)                              / q, std::exp(1.0)  },
        { "exp(-1)",pbf_exp_cf(-PBF_Q)                             / q, 1.0 / std::exp(1.0) },
        { "exp(2)", pbf_exp_cf(2 * PBF_Q)                          / q, std::exp(2.0)  },
    };

    std::printf("\n  %-10s %12s %12s %10s\n", "Function", "CF Result", "Exact", "Error");
    std::printf("  ------------------------------------------------\n");
    for (auto& t : tests) {
        double err = std::fabs(t.got - t.exact) / t.exact * 100.0;
        std::printf("  %-10s %12.6f %12.6f %9.4f%%\n", t.name, t.got, t.exact, err);
    }

    std::printf("\n=================================================================\n");
    std::printf("Pure CF. No tables. No floats in core ALU.\n");
    std::printf("=================================================================\n");
    return 0;
}
