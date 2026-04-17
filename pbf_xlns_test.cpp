// pbf_xlns_test.cpp — Sanity check that the xlns16 API shim works on top of PBF16.
// Runs the representative xlns16_* entrypoints and prints decoded results.

#include "pbf_xlns.cpp"

#include <cstdio>
#include <cmath>
#include <initializer_list>

static float relerr(float got, float expected) {
    if (fabsf(expected) < 1e-10f) return 0.0f;
    return fabsf(got - expected) / fabsf(expected) * 100.0f;
}

static int check(const char* label, float got, float expected, float tol_pct) {
    float err = relerr(got, expected);
    int ok = err <= tol_pct;
    std::printf("  %-30s got=%12.6f  expected=%12.6f  err=%6.3f%%  %s\n",
                label, got, expected, err, ok ? "OK" : "FAIL");
    return ok ? 0 : 1;
}

int main() {
    int failures = 0;
    std::printf("── xlns16 API shim over PBF16 ────────────────────────────────\n");

    // Conversion round-trip
    std::printf("\n[roundtrip]\n");
    for (float v : {1.0f, -1.0f, 2.0f, 0.5f, -0.25f, 1e3f, 1e-3f}) {
        char buf[32]; std::snprintf(buf, sizeof buf, "fp↔xlns16(%g)", v);
        failures += check(buf, xlns162fp(fp2xlns16(v)), v, 0.5f);
    }

    // Scalar arithmetic
    std::printf("\n[arithmetic]\n");
    xlns16 a = fp2xlns16( 3.0f);
    xlns16 b = fp2xlns16( 4.0f);
    // PBF16 mul/div SNR ≈ 33 dB (~2% error); add/sub ≈ 68 dB (~0.04%).
    failures += check("3 + 4",   xlns162fp(xlns16_add(a, b)), 7.0f, 1.0f);
    failures += check("3 - 4",   xlns162fp(xlns16_sub(a, b)), -1.0f, 5.0f);
    failures += check("3 * 4",   xlns162fp(xlns16_mul(a, b)), 12.0f, 3.0f);
    failures += check("12 / 3",  xlns162fp(xlns16_div(xlns16_mul(a, b), a)), 4.0f, 3.0f);
    failures += check("|-4|",    xlns162fp(xlns16_abs(xlns16_neg(b))), 4.0f, 0.5f);
    failures += check("sqrt(9)", xlns162fp(xlns16_sqrt(fp2xlns16(9.0f))), 3.0f, 0.5f);
    failures += check("3*4+1 (fma)", xlns162fp(xlns16_fma(a, b, xlns16_one)), 13.0f, 3.0f);

    // Constants
    std::printf("\n[constants]\n");
    failures += check("xlns16_one",     xlns162fp(xlns16_one),      1.0f,  0.5f);
    failures += check("xlns16_neg_one", xlns162fp(xlns16_neg_one), -1.0f,  0.5f);
    failures += check("xlns16_two",     xlns162fp(xlns16_two),      2.0f,  0.5f);
    failures += check("xlns16_half",    xlns162fp(xlns16_half),     0.5f,  0.5f);

    // Predicates
    std::printf("\n[predicates]\n");
    std::printf("  is_zero(0)     = %d (expect 1)\n", xlns16_is_zero(xlns16_zero));
    std::printf("  is_negative(-1)= %d (expect 1)\n", xlns16_is_negative(xlns16_neg_one));
    std::printf("  is_positive(2) = %d (expect 1)\n", xlns16_is_positive(xlns16_two));

    // Comparisons
    std::printf("\n[comparisons]\n");
    std::printf("  3 < 4  : %d (expect 1)\n", xlns16_lt(a, b));
    std::printf("  3 > 4  : %d (expect 0)\n", xlns16_gt(a, b));
    std::printf("  3 == 3 : %d (expect 1)\n", xlns16_eq(a, a));

    // Transcendentals
    std::printf("\n[transcendentals]\n");
    failures += check("exp(1)",  xlns162fp(xlns16_exp(xlns16_one)), std::exp(1.0f), 1.0f);
    failures += check("log(e)",  xlns162fp(xlns16_log(fp2xlns16(std::exp(1.0f)))), 1.0f, 1.0f);
    failures += check("2^3",     xlns162fp(xlns16_pow(xlns16_two, fp2xlns16(3.0f))), 8.0f, 1.0f);

    // Activations
    std::printf("\n[activations]\n");
    failures += check("relu(-1)",  xlns162fp(xlns16_relu(xlns16_neg_one)),  0.0f, 0.0f);
    failures += check("relu(+2)",  xlns162fp(xlns16_relu(xlns16_two)),       2.0f, 0.5f);
    failures += check("tanh(0)",   xlns162fp(xlns16_tanh(xlns16_zero)),      0.0f, 0.0f);
    failures += check("sigmoid(0)",xlns162fp(xlns16_sigmoid(xlns16_zero)),   0.5f, 0.5f);

    // Vector ops
    std::printf("\n[vector]\n");
    xlns16 va[4] = { fp2xlns16(1.0f), fp2xlns16(2.0f), fp2xlns16(3.0f), fp2xlns16(4.0f) };
    xlns16 vb[4] = { fp2xlns16(1.0f), fp2xlns16(1.0f), fp2xlns16(1.0f), fp2xlns16(1.0f) };
    failures += check("sum([1..4])", xlns162fp(xlns16_sum(va, 4)), 10.0f, 1.0f);
    failures += check("dot([1..4],1)", xlns162fp(xlns16_vec_dot(va, vb, 4)), 10.0f, 3.0f);

    // Batch ops
    std::printf("\n[batch]\n");
    xlns16 vc[4];
    xlns16_batch_mul(va, vb, vc, 4);
    float check_sum = 0.0f;
    for (int i = 0; i < 4; i++) check_sum += xlns162fp(vc[i]);
    failures += check("Σ batch_mul", check_sum, 10.0f, 3.0f);

    std::printf("\n══════════════════════════════════════════════════════════════\n");
    std::printf("Result: %s  (%d failure%s)\n",
                failures == 0 ? "PASS" : "FAIL",
                failures, failures == 1 ? "" : "s");
    return failures == 0 ? 0 : 1;
}
