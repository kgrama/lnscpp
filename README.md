# pbf16

Table-free 16-bit LNS arithmetic for C++. Drop-in replacement for [`xlns16.cpp`](https://github.com/xlnsresearch/xlnscpp) — same types, same function names, same operator overloads.

*Patent pending GB 2602876.1*

---

## Files

| File | Description |
|------|-------------|
| `pbf16.cpp` | PBF encoding — CF-based sb/db, no lookup tables |
| `pbf16_sbp.cpp` | SBP encoding — NOT(x) = −x, CLZ-adaptive CF depth |
| `tests/sbp_vs_pbf_test.cpp` | Comparison test: PBF vs SBP vs bfloat16 |
| `Makefile` | Builds and runs all three |

---

## Usage  -- Broken

Replace:
```cpp
#define xlns16_ideal
#include "xlns16.cpp"
```
With:
```cpp
#include "pbf16.cpp"   // or pbf16_sbp.cpp
```

Remove `xlns16cvtbl.h`, `xlns16sbdbtbl.h`, `xlns16logtbl.h`, `xlns16exptbl.h`.

Requires C++11 with `__int128` (gcc/clang, x86-64 or ARM64).

---

## API

```cpp
xlns16  fp2xlns16(float x)
float   xlns162fp(xlns16 x)
xlns16  xlns16_add(xlns16 a, xlns16 b)
xlns16  xlns16_sub(xlns16 a, xlns16 b)
xlns16  xlns16_mul(xlns16 a, xlns16 b)
xlns16  xlns16_div(xlns16 a, xlns16 b)
xlns16  xlns16_neg(xlns16 a)
xlns16  xlns16_abs(xlns16 a)
```

Vector, batch, and activation functions (`xlns16_vec_dot`, `xlns16_sigmoid`, `xlns16_softmax`, etc.) also present.

---

## Encoding

**PBF** (`pbf16.cpp`): sign bit at position 15 (XOR mask). Identical internal format to xlns16.

**SBP** (`pbf16_sbp.cpp`): Sign-By-Position — values below mirror axis are negative, above are positive. Negation is bitwise NOT (single gate): `~x == neg(x)`.

```
SBP layout (16-bit):
  0x0000        = Zero
  0x0001–0x7FFF = Negative values
  0x8000–0xFFFE = Positive values
  0xFFFF        = NaN  (= ~Zero)
```

---

## How sb/db are computed

LNS addition requires the Gauss log functions:

- `sb(z) = log₂(1 + 2⁻ᶻ)`
- `db(z) = log₂(1 − 2⁻ᶻ)`

Both files evaluate these via continued fraction over Q32 fixed-point integers using `__int128`. CF depth is CLZ-adaptive: deeper evaluation when operands are close in magnitude.

A 16-bit ROM table requires ~32K×16 bits of storage (~524 KB), which is impractical in gate-array or table-free silicon contexts. The CF approach trades ROM for combinational logic; no gate-count synthesis has been performed in this repo.

---

## Measured SNR

Random inputs, range ±16, N=20000:

| Format     | add    | sub    | mul    | div    | round-trip | cancellation |
|------------|--------|--------|--------|--------|------------|--------------|
| PBF        | 53.0   | 53.1   | 53.0   | 52.2   | 76.5       | 41.0 dB      |
| SBP        | 53.0   | 53.1   | 53.0   | 52.2   | 76.5       | 41.0 dB      |
| bfloat16   | 46.3   | 46.4   | 41.7   | 48.3   | 67.0       | 39.7 dB      |

PBF and SBP share the same CF engine; SNR is identical. bfloat16 cancellation test: `a + (−a×(1+ε))`, ε ∈ [0.001, 0.5].

---

## Running the tests

```bash
make        # builds pbf_test, sbp_test, bfloat_test and runs all three
make clean
```

---

## Differences from xlns16

| Feature              | xlns16          | pbf16 / pbf16_sbp     |
|----------------------|-----------------|-----------------------|
| sb/db method         | Table or float  | CF, Q32 fixed-point   |
| ROM required         | Yes (5–4000 KB) | None                  |
| `__int128` required  | No              | Yes                   |
| Compile-time options | ideal/alt/table | None (CF always)      |
| NOT(x) = −x          | No              | pbf16_sbp only        |

---

## License

Research and evaluation use only. Commercial use prohibited without a separate written license.

Patent pending GB 2602876.1. No patent license is granted by the software license.

See [LICENSE](LICENSE) for full terms.
