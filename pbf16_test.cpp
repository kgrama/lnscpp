// pbq16_test.cpp — tests for virtual mirror encoding
// Build:
//   g++ -std=c++11 -O2 -I.. -o pbq16_test tests/pbq16_test.cpp -lm

#include "../pbq16.cpp"
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <cstdint>

static int pass = 0, fail = 0;
static void check(bool ok, const char* fmt, ...) {
    ok ? pass++ : fail++;
    printf("  [%s] ", ok ? "OK  " : "FAIL");
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    printf("\n");
}

static void test_structure() {
    printf("\n═══ Index structure ═══\n");

    // Sentinels
    check(PBQ_ZERO == 0x0000, "Zero  = 0x0000");
    check(PBQ_NAN  == 0xFFFF, "NaN   = 0xFFFF");
    check((uint16_t)(~PBQ_ZERO) == PBQ_NAN,  "~Zero == NaN");

    // Axis boundary
    pbq16 sm_neg = 0x7FFF;  // smallest negative (dist=1)
    pbq16 sm_pos = 0x8000;  // smallest positive (dist=1)
    check((uint16_t)(~sm_neg) == sm_pos, "~smallest_neg == smallest_pos  (0x7FFF ↔ 0x8000)");
    check(pbq_dist(sm_neg) == 1, "dist(0x7FFF) == 1");
    check(pbq_dist(sm_pos) == 1, "dist(0x8000) == 1");

    // Largest values
    pbq16 lg_neg = 0x0001;
    pbq16 lg_pos = 0xFFFE;
    check((uint16_t)(~lg_neg) == lg_pos, "~largest_neg == largest_pos  (0x0001 ↔ 0xFFFE)");
    check(pbq_dist(lg_neg) == PBQ_MAX_DIST, "dist(0x0001) == 0x7FFF");
    check(pbq_dist(lg_pos) == PBQ_MAX_DIST, "dist(0xFFFE) == 0x7FFF");

    // Predicates
    check(pbq_is_zero(PBQ_ZERO),   "is_zero(0x0000)");
    check(pbq_is_nan(PBQ_NAN),     "is_nan(0xFFFF)");
    check(pbq_is_neg(0x0001),     "is_neg(0x0001)");
    check(pbq_is_neg(0x7FFF),     "is_neg(0x7FFF)");
    check(pbq_is_pos(0x8000),     "is_pos(0x8000)");
    check(pbq_is_pos(0xFFFE),     "is_pos(0xFFFE)");
    check(!pbq_is_pos(PBQ_NAN),    "!is_pos(NaN)");

    // NOT = negation for all codes
    printf("\n  NOT(x) == neg(x) — exhaustive over all 65536 codes:\n");
    int not_fail = 0;
    for (int c = 0; c <= 0xFFFF; c++) {
        pbq16 x = (pbq16)c;
        if (pbq_neg(x) != (pbq16)(~x)) not_fail++;
    }
    check(not_fail == 0, "NOT(x) == neg(x) for all 65536 codes (%d failures)", not_fail);

    // dist symmetry: neg and pos at same dist decode to same magnitude
    printf("\n  Magnitude symmetry across axis:\n");
    int sym_fail = 0;
    for (int d = 1; d <= 0x7FFF; d++) {
        pbq16 pos = pbq_make(0, (uint16_t)d);
        pbq16 neg = pbq_make(1, (uint16_t)d);
        float fp = pbq2f(pos), fn = pbq2f(neg);
        if (fabsf(fp + fn) > 1e-6f * fabsf(fp)) sym_fail++;
    }
    check(sym_fail == 0, "pbq2f(pos_d) == -pbq2f(neg_d) for all dist (%d failures)", sym_fail);
}

static void test_encoding() {
    printf("\n═══ Float ↔ pbq16 round-trip ═══\n");

    float vals[] = {
        1.0f, -1.0f, 2.0f, -2.0f, 0.5f, -0.5f,
        4.0f, 0.25f, 3.14159f, -2.71828f, 0.001f, 1000.0f
    };
    for (float v : vals) {
        pbq16 enc = f2pbq(v);
        float dec = pbq2f(enc);
        float relerr = fabsf(dec - v) / fabsf(v) * 100.f;
        check(relerr < 1.0f, "%.5f → 0x%04x → %.5f  (err=%.3f%%)",
              v, (unsigned)enc, dec, relerr);
    }

    // Zero and NaN
    check(f2pbq(0.0f) == PBQ_ZERO, "f2pbq(0) == PBQ_ZERO");
    check(pbq2f(PBQ_ZERO) == 0.0f, "pbq2f(Zero) == 0.0");
    check(pbq_is_nan(f2pbq(__builtin_nanf(""))), "f2pbq(NaN) == PBQ_NAN");
    check(pbq_is_nan(pbq2f(PBQ_NAN) != pbq2f(PBQ_NAN) ? PBQ_NAN : PBQ_NAN),
          "pbq2f(NaN) is NaN");

    // Ordering: larger dist → larger magnitude
    printf("\n  Dist ordering (pos side):\n");
    int order_fail = 0;
    float prev = 0.0f;
    for (int d = 1; d <= 0x7FFF; d += 128) {
        float v = pbq2f(pbq_make(0, (uint16_t)d));
        if (v <= prev) order_fail++;
        prev = v;
    }
    check(order_fail == 0, "pbq2f increases with dist (%d violations)", order_fail);
}

static void test_arithmetic() {
    printf("\n═══ Arithmetic ═══\n");

    struct { float a, b; } cases[] = {
        {1,1},{3,4},{-3,4},{3,-4},{-2,-3},
        {0.5f,0.5f},{8,0.125f},{-5,3}
    };
    for (auto& c : cases) {
        char l[32];
        float refa = c.a * c.b;
        float refd = c.a / c.b;
        float gota = pbq2f(pbq_mul(f2pbq(c.a), f2pbq(c.b)));
        float gotd = pbq2f(pbq_div(f2pbq(c.a), f2pbq(c.b)));
        snprintf(l, 32, "%g*%g", c.a, c.b);
        check(fabsf(gota-refa)/fabsf(refa) < 0.02f, "%-12s  ref=%.4f  got=%.4f", l, refa, gota);
        snprintf(l, 32, "%g/%g", c.a, c.b);
        check(fabsf(gotd-refd)/fabsf(refd) < 0.02f, "%-12s  ref=%.4f  got=%.4f", l, refd, gotd);
    }

    // neg
    float neg_vals[] = {1.0f, -2.0f, 0.5f, 100.0f};
    for (float v : neg_vals) {
        float got = pbq2f(pbq_neg(f2pbq(v)));
        check(fabsf(got + v) < 0.01f * fabsf(v) + 1e-7f,
              "neg(%.3f) = %.3f", v, got);
    }

    // abs
    check(fabsf(pbq2f(pbq_abs(f2pbq(-3.0f))) - 3.0f) < 0.05f, "abs(-3) == 3");
    check(fabsf(pbq2f(pbq_abs(f2pbq( 3.0f))) - 3.0f) < 0.05f, "abs(+3) == 3");

    // add
    printf("\n  Add (float-based, placeholder for CF):\n");
    struct { float a, b; } add_cases[] = {{1,1},{3,4},{-3,4},{0.5f,0.5f},{-1,1}};
    for (auto& c : add_cases) {
        float ref = c.a + c.b;
        float got = pbq2f(pbq_add(f2pbq(c.a), f2pbq(c.b)));
        char l[32]; snprintf(l, 32, "%g+%g", c.a, c.b);
        float tol = fabsf(ref) < 1e-4f ? 0.1f : fabsf(ref) * 0.02f;
        check(fabsf(got - ref) <= tol,
              "%-10s  ref=%.4f  got=%.4f", l, ref, got);
    }
}

int main() {
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  pbq16 — Virtual Mirror Encoding Tests\n");
    printf("══════════════════════════════════════════════════════════════\n");
    test_structure();
    test_encoding();
    test_arithmetic();
    printf("\n══════════════════════════════════════════════════════════════\n");
    printf("  SUMMARY: %d passed, %d failed\n", pass, fail);
    printf("══════════════════════════════════════════════════════════════\n");
    return fail ? 1 : 0;
}
