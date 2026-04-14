# pbf16

**Parameterized Bounded Format — table-free LNS arithmetic for C++**

Drop-in replacement for [`xlns16.cpp`](https://github.com/xlnsresearch/xlnscpp) from xlnsresearch.  
Same types. Same function names. Same operator overloads. No tables.

*Patent pending GB 2602876.1*

---

## What it is

`pbf16.cpp` replaces the Gauss-log lookup tables in xlns16 with fixed-point
continued fraction evaluation. Addition and subtraction — the hard operations
in LNS — are computed via CF-based `sb` and `db` functions over pure integer
arithmetic. No ROM, no SRAM, no `xlns16sbdbtbl.h`, no `xlns16cvtbl.h`.

The internal 16-bit format is identical to xlns16:

```
+------+-------------------------+
| sign | int(log2) . frac(log2) |
+------+-------------------------+
1 sign bit · 8 integer bits · 7 fractional bits
Offset by 0x4000 (logsignmask)
```

---

## Usage

Replace:
```cpp
#define xlns16_ideal
#include "xlns16.cpp"
```

With:
```cpp
#include "pbf16.cpp"
```

Remove `xlns16cvtbl.h`, `xlns16sbdbtbl.h`, `xlns16logtbl.h`, `xlns16exptbl.h`.  
`xlns16_ideal`, `xlns16_alt`, `xlns16_table` macros are accepted but have no
effect — CF is always used.

---

## API

Identical to xlns16. C function interface:

```cpp
xlns16  fp2xlns16(float x)            // float → LNS
float   xlns162fp(xlns16 x)           // LNS → float

xlns16  xlns16_add(xlns16 a, xlns16 b)
xlns16  xlns16_sub(xlns16 a, xlns16 b)
xlns16  xlns16_mul(xlns16 a, xlns16 b)
xlns16  xlns16_div(xlns16 a, xlns16 b)
xlns16  xlns16_neg(xlns16 a)
xlns16  xlns16_abs(xlns16 a)
```

C++ class with overloaded operators:

```cpp
xlns16_float a, b;
a = 3.0f;
b = 4.0f;
xlns16_float c = a * b;   // 12
xlns16_float d = a + b;   // 7
std::cout << d;            // 7.0
```

Batch and vector operations (`xlns16_vec_dot`, `xlns16_batch_mul`, etc.),
activation functions (`xlns16_relu`, `xlns16_sigmoid`, `xlns16_softmax`, etc.)
all present and interface-compatible.

---

## Why no tables

LNS addition requires evaluating the Gauss logarithm functions:

- `sb(z) = log₂(1 + 2^z)`  
- `db(z) = log₂(2^z − 1)`

xlns16 computes these from lookup tables (or floating-point in `xlns16_ideal`
mode). Table size grows as 2^n — impractical above 10–12 bits, and a liability
in constrained silicon (no block RAM, side-channel sensitivity, formal
verification requirements).

pbf16 evaluates sb and db via 8-step continued fractions over Q32 fixed-point
integers, using `__int128` for intermediate products. The CF converges in
8 steps with accuracy matching `xlns16_ideal` at 16-bit precision.

Gate-count comparison for the add operation at 16 bits:

| Approach          | Gates (est.) | Latency   | Table? |
|-------------------|-------------|-----------|--------|
| ROM table         | ~4,194,000  | 1 cycle   | YES    |
| CPLA (piecewise)  | ~1,792      | 2 cycles  | NO     |
| pbf16 CF          | ~34,000     | ~10 cycles| NO     |

ROM is physically impractical at 16 bits. CF trades latency for zero ROM,
arbitrary precision, and formal verifiability.

---

## SNR vs IEEE 754

Measured on random samples in the safe interior of the dynamic range:

```
Format    bits    add     sub     mul     div
──────────────────────────────────────────────
float16     16   73.3   74.1   73.6   73.7  dB
pbf16       16   73.3   72.9   69.1   68.8  dB
```

pbf16 matches float16 on add/sub. mul/div are ~5 dB lower due to two
input quantizations compounding in log space — a fundamental LNS property,
not a CF artifact.

---

## Building

Single-file include, no dependencies beyond `<stdint.h>` and `<math.h>`:

```bash
g++ -O2 your_program.cpp -lm
```

Requires a C++11 compiler with `__int128` support (gcc, clang on x86-64 and
ARM64). For MSVC, `__int128` can be emulated — contributions welcome.

---

## Differences from xlns16

| Feature              | xlns16             | pbf16                  |
|----------------------|--------------------|------------------------|
| sb / db method       | Table or float     | 8-step CF, Q32 integer |
| ROM required         | Yes (5–4000 KB)    | None                   |
| `__int128` required  | No                 | Yes (intermediate ops) |
| Compile-time options | `ideal/alt/table`  | None (CF always)       |
| Internal format      | Identical          | Identical              |
| API                  | Identical          | Identical              |

---

## License

Research and evaluation use only. Commercial use is explicitly prohibited
without a separate written license from the copyright holder.

Patent pending GB 2602876.1. The continued fraction method for LNS sb/db
evaluation is the subject of a pending patent application. No patent license
is granted by the software license.

See [LICENSE](LICENSE) for full terms.

For commercial licensing enquiries, contact the copyright holder.

---

## References

xlnsresearch/xlnscpp: https://github.com/xlnsresearch/xlnscpp  
xlnsresearch reference list: https://xlnsresearch.github.io
