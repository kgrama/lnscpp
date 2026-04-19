# lnscpp

An exact, table-free logarithmic number system (LNS) in C++. The format
is **PBF** (Parameterized Bounded Format): one runtime-parameterised
encoding that covers PBF4 through PBF16 and arbitrary custom widths
through `pbf_init(&p, n_bits, v_min, v_max)`.

## Why this is different

Two properties separate PBF from traditional LNS implementations
(`xlns16`, bfloat16-style) and from float:

- **Exact `sb` / `db`.** The Gauss log functions
  `sb(d) = log(1 + e^d)` and `db(d) = log|1 − e^d|` are evaluated by a
  closed-form continued fraction in Q32 fixed-point. Results are
  correctly rounded at the working precision — no ROM, no interpolation
  table, no error budget beyond Q32 rounding. The same Q32 engine
  serves every bit width, so format choice is a runtime parameter
  rather than a separate implementation per width.

- **Sign-By-Position (SBP) encoding.** No sign bit. Sign is determined
  by where the code sits relative to the midpoint, so negation is
  bitwise complement (`max_code − a`) and unsigned integer comparison
  matches value comparison. Zero and infinity are the fixed points;
  ±0 ambiguity does not arise. The encoding and its structural proofs
  are in [docs/SBP_ENCODING.md](docs/SBP_ENCODING.md) and
  [proofs/](proofs/).

The trade-off vs float is the usual LNS one: mul/div become an integer
add, add/sub become a fixed-depth CF evaluation, and per-bit precision
is lower (PBF8 ≈ 4.6 bits, PBF16 ≈ 11.7 bits).

---

## Quickstart

```cpp
// hello.cpp — compile:  g++ -std=c++11 -O2 hello.cpp -lm
#include "pbf.cpp"
#include <cstdio>

int main() {
    pbf_t p;
    pbf16_init(&p);                          // 16-bit, v_min=1e-8, v_max=1e8

    uint32_t a = pbf_encode(&p, 3.14);
    uint32_t b = pbf_encode(&p, 2.0);
    uint32_t c = pbf_mul(&p, a, b);          // integer add on log levels
    uint32_t d = pbf_add(&p, c, pbf_encode(&p, 1.0));

    printf("3.14 * 2.0 + 1.0 = %f\n", pbf_decode(&p, d));
    return 0;
}
```

## Presets

| Preset             | Bits | Range          | Effective precision           | Typical use                                  |
|--------------------|------|----------------|-------------------------------|----------------------------------------------|
| `pbf4_init`        | 4    | 0.1 … 10       | ~2 bits                       | demos, didactic examples                     |
| `pbf8_init`        | 8    | 10⁻⁴ … 10⁴     | ~4.6 bits                     | aggressive NN weight/activation quantization |
| `pbf16_init`       | 16   | 10⁻⁸ … 10⁸     | ~11.7 bits                    | drop-in for bfloat16-level precision work    |
| `pbf_taper8_init`  | 8    | 10⁻⁴ … 10⁴     | ~5.4 bits at unity            | tapered: dense around \|x\|=1.0              |
| `pbf_taper16_init` | 16   | 10⁻⁸ … 10⁸     | ~13.8 bits at unity, ~11 wide | tapered: dense around \|x\|=1.0              |

Custom ranges via `pbf_init(&p, n_bits, v_min, v_max)` or
`pbf_taper_init(&p, n_bits, v_max)` (tapered uses `[1/v_max, v_max]`).

---

## Requirements

- C++11 compiler with `__int128` (gcc or clang on x86-64 / ARM64,
  Linux / macOS / WSL). MSVC does **not** support `__int128` — on
  Windows use gcc under WSL or MSYS2.
- `make`, `libm`.

## Build

```bash
make          # build everything (PBF + xlns32)
make pbf      # PBF only
make xlns32   # xlns32 only
make run      # build and run the PBF demo + shim test
make clean    # remove build/
```

All binaries land in `build/`. To run a single test:

```bash
build/pbftest           # SNR + CF-accuracy demo
build/pbf_taper_test    # tapered vs flat SNR, split by operand regime
build/pbf_bench         # speed bench: PBF vs tapered PBF vs xlns32
build/pbf_xlns_test     # xlns16 API shim verification
```

---

## API

### Scalar (from `pbf.cpp`)

```cpp
void     pbf_init   (pbf_t* p, int n_bits, double v_min, double v_max);
uint32_t pbf_encode (const pbf_t* p, double v);
double   pbf_decode (const pbf_t* p, uint32_t code);
uint32_t pbf_neg    (const pbf_t* p, uint32_t a);
uint32_t pbf_add    (const pbf_t* p, uint32_t a, uint32_t b);
uint32_t pbf_sub    (const pbf_t* p, uint32_t a, uint32_t b);
uint32_t pbf_mul    (const pbf_t* p, uint32_t a, uint32_t b);
uint32_t pbf_div    (const pbf_t* p, uint32_t a, uint32_t b);
int      pbf_similarity(const pbf_t* p, uint32_t a, uint32_t b);
```

### Batch (from `pbf_batch.cpp`)

```cpp
void pbf_batch_add   (const pbf_t*, const uint32_t* a, const uint32_t* b, uint32_t* c, size_t n);
// … _sub, _mul, _div, _scale, _neg, _abs, _encode, _decode

uint32_t pbf_sum     (const pbf_t*, const uint32_t* a, size_t n);
uint32_t pbf_dot     (const pbf_t*, const uint32_t* a, const uint32_t* b, size_t n);
uint32_t pbf_max     (const pbf_t*, const uint32_t* a, size_t n);
uint32_t pbf_min     (const pbf_t*, const uint32_t* a, size_t n);
double   pbf_vec_dot_f64(const pbf_t*, const double* a, const double* b, size_t n);
```

### Tapered (from `pbf_taper.cpp`)

```cpp
void     pbf_taper_init  (pbf_taper_t* p, int n_bits, double v_max);
uint32_t pbf_taper_encode(const pbf_taper_t* p, double v);
double   pbf_taper_decode(const pbf_taper_t* p, uint32_t code);
uint32_t pbf_taper_neg   (const pbf_taper_t* p, uint32_t a);
uint32_t pbf_taper_add   (const pbf_taper_t* p, uint32_t a, uint32_t b);
uint32_t pbf_taper_sub   (const pbf_taper_t* p, uint32_t a, uint32_t b);
uint32_t pbf_taper_mul   (const pbf_taper_t* p, uint32_t a, uint32_t b);
uint32_t pbf_taper_div   (const pbf_taper_t* p, uint32_t a, uint32_t b);
```

Four tapers with step ratios 8:1:1:8 — inner tapers (around |x|=1.0)
have 8× finer log step than the outer tapers. The W16 boundaries
(8191, 16382, 24573, 32766) match the partition proven in
[proofs/P32_Taper.v](proofs/P32_Taper.v); the per-taper step sizes are
an interpretation-layer choice and don't affect the structural soundness
proof. Arithmetic round-trips through Q32 log space, so mul/div pay an
encode/decode pass vs flat PBF (see Speed bench below).

### xlns16 shim (from `pbf_xlns.cpp`)

Drop-in for code written against the `xlns16_*` naming convention
(`fp2xlns16`, `xlns16_add`, `xlns16_one`, activations, softmax,
layernorm, …), backed by PBF16. A lazily-initialised `pbf_t` context
sits behind the API. Decoded values match; bit patterns do not.

---

## SNR

Flat PBF, full-range operand sweep (`./build/pbftest`):

```
Format       Roundtrip    Add        Sub        Mul        Div
PBF8         28.0 dB      26.1 dB    27.8 dB    20.7 dB    20.4 dB
PBF16        70.7 dB      68.0 dB    69.8 dB    64.7 dB    63.9 dB
```

Tapered PBF16, per operand regime (`./build/pbf_taper_test`):

```
Regime                    Roundtrip   Mul       Div       Add
near-unity [0.5, 2]        82.9 dB    79.9 dB   79.6 dB   57.1 dB
mid-range  [1e-2, 1e2]     64.9 dB    59.4 dB   59.8 dB   60.9 dB
wide       [1e-7, 1e7]     64.6 dB    58.2 dB   59.1 dB   64.1 dB
```

Tapered buys +13 dB roundtrip / +17 dB mul near unity vs flat PBF16,
at the cost of ~6 dB at the extremes. Useful when operands cluster
around 1.0 (NN activations, normalised signals).

Rule of thumb: ~6 dB per bit. Flat PBF mul/div trail add/sub by ~5 dB
because the result composes three independent rounding operations
(`encode(x)`, `encode(y)`, and the integer offset) where add/sub
composes only one final round of the CF output.

The CF `ln`, `exp`, `sb`, `db` primitives are correctly rounded at Q32;
`./build/pbftest` prints agreement with `libm` to its full displayed
precision.

---

## Speed

Single-thread throughput, pre-encoded operands (`./build/pbf_bench`):

```
                       add       sub       mul       div
PBF16               94.5 ns   94.7 ns    2.5 ns    2.7 ns
Tapered PBF16      105.3 ns  106.4 ns   13.1 ns   12.8 ns
xlns32 (LUT)         9.2 ns    9.4 ns    2.2 ns    2.2 ns
```

- **Mul/div: parity with xlns32.** All three formats reduce mul/div to
  an integer add on log levels; the 10–20 % spread is bounds-checks
  and the offset add.
- **Add/sub: xlns32 wins ~10×.** xlns32 looks up `sb`/`db` from a
  ~12 KB ROM (3-stage table); PBF runs the Q32 continued fraction
  (~30 muldivs per call). PBF buys deterministic timing, no LUT, no
  cache pollution, and runtime-parameterised bit width — at 10× the
  cycle count.
- **Tapered mul/div: ~5× the flat-PBF cost.** Each operand round-trips
  through Q32 log space (decode → integer add → encode); add/sub barely
  move because the CF already dominates. A 64 KB precomputed
  `mag → log` LUT would close most of the mul/div gap, at the cost of
  the table-free property.

---

## Fractal demo

![Log Mirror fractal, PBF-encoded complex iteration](log-mirror.png)

Open [`fractal_zoo.html`](fractal_zoo.html) directly in a browser — no
server needed. Requires WebGL 1 (Chrome 80+, Firefox 51+, Safari 15+,
Edge, mobile browsers).

Seven fractals (Mandelbrot, Julia, Burning Ship, Tricorn, Log Mirror,
Lyapunov, Newton z³−1) iterate through an integer PBF ALU implemented
in GLSL ES 1.0:

- complex values are carried as `ivec2 = (re_code, im_code)`
- `pbf_mul` / `pbf_div` are integer adds on magnitude levels
- `pbf_add` uses `log(1+exp(d))` at the float boundary (WebGL 1 lacks
  int64 for the Q32 CF path used in `pbf.cpp`)
- Lyapunov reads `log|deriv|` directly off the PBF level — no `log()`
  in the inner loop

A sidebar toggle switches between **PBF8** (chunky; ~4.6 bits) and
**PBF16** (smooth; ~11.7 bits). Newton and Log Mirror chain enough
divisions per iteration that PBF8 typically fails to converge; PBF16
works for all fractals.

---

## Format details

A PBF code is an unsigned integer in `[0, 2ⁿ − 1]`:

```
code == 0            → zero
code <  mid          → negative (magnitude grows as code decreases)
code >= mid          → positive (magnitude grows as code increases)
code == max_code     → infinity
```

Negation is `pbf_neg(a) = max_code − a` (ones' complement within the
range); ordering is unsigned integer ordering. The `pbf_t` descriptor
holds `n_bits`, `v_min`, `v_max`, and precomputed `scale` / `offset`.

Internally the ALU is Q32 fixed-point. `ln`, `exp`, `sb`, `db` are
evaluated with continued fractions; `__int128` is used only for the
intermediate multiplies in `muldiv`.

The encoding's structural claims — mirror symmetry, negation by
complement, partition totality, ordering, minimum-magnitude bound — are
written up in [docs/SBP_ENCODING.md](docs/SBP_ENCODING.md) and
mechanically verified in [proofs/](proofs/) (Coq/Rocq).

---

## Files

### PBF

| File                | Contents                                              |
|---------------------|-------------------------------------------------------|
| `pbf.cpp`            | Core: Q32 CF engine, encode/decode, neg/add/sub/mul/div |
| `pbf_taper.cpp`      | Tapered PBF (4 tapers, denser around \|x\|=1.0)        |
| `pbf_batch.cpp`      | Batch, vector, and reduction ops over PBF arrays       |
| `pbf_xlns.cpp`       | `xlns16_*` API shim backed by PBF16                    |
| `pbftest.cpp`        | SNR demo, error distribution, CF-primitive accuracy    |
| `pbf_taper_test.cpp` | Tapered-vs-flat SNR by operand regime                  |
| `pbf_bench.cpp`      | Speed bench: PBF / Tapered PBF / xlns32                |
| `pbf_xlns_test.cpp`  | Shim verification                                      |
| `pbf.py`             | Python reference used for the C++ port                 |
| `fractal_zoo.html`   | WebGL demo — every fractal runs the PBF integer ALU in the shader |

### Vendored: xlns32

32-bit LNS from the `xlnsresearch/xlnscpp` project. Included unchanged
as a reference and test target; PBF does not depend on it.

| File                                 | Contents                  |
|--------------------------------------|---------------------------|
| `xlns32.cpp`, `xlns32tbl.h`          | 32-bit LNS implementation |
| `xlns32test.cpp`                     | Numeric test suite        |
| `xlns32funtest.cpp`                  | Transcendental functions  |
| `xlns32_new_functions_test.cpp`      | Batch/activation tests    |
| `tests/xlns32_*.cpp`                 | Focused feature tests     |

---

## License

Research and evaluation use only. Commercial use prohibited without a
separate written license. Patent pending GB 2602876.1. No patent license
is granted by the software license. See [LICENSE](LICENSE) for full terms.
