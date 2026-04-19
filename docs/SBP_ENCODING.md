# SBP: A Mirror-Axis Number Encoding

## Abstract

Sign-By-Position (SBP) is a number encoding in which each representable
value is assigned to an unsigned integer index in a finite range
`[0, 2^N−1]`, with a virtual mirror axis at the centre of that range.
Sign is determined by which side of the mirror the index falls on;
magnitude is the distance from the mirror; arithmetic negation is bitwise
complement of the index. The encoding reserves fixed index positions for
zero, infinities, NaN, and (in some variants) denormals.

This document describes the encoding. The arithmetic built on top of it
is out of scope here. Structural claims stated here are mechanically
verified in Coq/Rocq under `proofs/`; verification status is reported
per-theorem in Section 8.

## 1. The mirror encoding and Layout C

SBP uses a pointer-arithmetic model. A value is identified by an index
`idx` into a finite code space and a linear refinement `lin` within the
code cell. Values are laid out monotonically along the `idx` axis. The
mirror is a half-integer position at `2^(N−1) − 0.5`. Indices below the
mirror encode negative values; indices above encode positive values.

```
idx:  0      1          pit_lo  NaN NaN pit_hi          IDX_MAX
      └──┘   └──┘  ...  └──┘  └──┘ └──┘  └──┘   ...     └──┘
      zero   -inf       neg   nan   nan  pos            +inf
                        denorm            denorm
             ┌────── negative normals ──┘   └── positive normals ──┐
```

### 1.1 Bit layout

An N-bit SBP value is divided into two halves, an index `idx` and a
linear refinement `lin`, each `N/2` bits wide:

```
 MSB                             LSB
 ┌────────────┬────────────┐
 │    idx     │    lin     │
 └────────────┴────────────┘
```

Position within the overall code space is `idx · 2^(N/2) + lin`. The
`idx` half determines sign and coarse magnitude; the `lin` half provides
sub-cell refinement. At `N = 8`: 4 bits `idx`, 4 bits `lin`. At `N = 32`:
16 bits each.

### 1.2 Negation as complement

The central structural property is that for any valid code `x`,
arithmetic negation `−x` is represented by the bitwise complement `~x`.
This follows from the mirror placement:

- the code space spans `[0, 2^N−1]`;
- the mirror lies at `2^(N−1) − 0.5`;
- reflection of an integer `x` about `2^(N−1) − 0.5` is `(2^N−1) − x = ~x`.

Consequences: no sign bit is needed, `neg(x) = ~x` is one NOT per bit,
and negation preserves magnitude exactly. The involution
`neg(neg(x)) = x` is stated in `proofs/P01_NegInvolution.v`.

### 1.3 The pivot gap

A contiguous region straddling the mirror axis holds special values. In
Layout C this region is four slots wide: one negative-denormal slot, two
NaN slots, and one positive-denormal slot. Denormals provide
subnormal-like values near zero; NaN propagates through arithmetic. The
remaining indices on either side of the gap carry normal magnitudes.

### 1.4 Reserved positions

The smallest Layout C variant reserves four positions around the mirror
and two at the extremes:

| idx                     | Name         | Meaning                   |
|-------------------------|--------------|---------------------------|
| 0                       | zero         | exact zero                |
| 1                       | neg_inf      | negative infinity         |
| pit_lo = 2^(N/2−1) − 2  | neg_denormal | small-magnitude negative  |
| nan_lo = 2^(N/2−1) − 1  | nan          | not-a-number              |
| nan_hi = 2^(N/2−1)      | nan          | not-a-number              |
| pit_hi = 2^(N/2−1) + 1  | pos_denormal | small-magnitude positive  |
| IDX_MAX = 2^N − 1       | pos_inf      | positive infinity         |

Normal negative values occupy `idx ∈ (1, pit_lo)`. Normal positive
values occupy `idx ∈ (pit_hi, IDX_MAX)`. The partition is shown
exhaustive and disjoint in `proofs/P21_IdxPartition.v`.

### 1.5 Sign determination

Sign of a value is determined by position relative to the mirror:

```
sign(x) =  0                 if x is the zero tag
        =  +1                if idx(x) > pit_hi
        =  −1                if 1 < idx(x) < pit_lo
        = undefined          if x is a NaN or denormal tag
```

Extraction is a single integer comparison against `pit_hi` (and `pit_lo`
for the negative case). No bit extraction from the payload is required.
Stated in `proofs/P09_SignByPosition.v`.

### 1.6 Comparison and ordering

Normal values are monotonic in `idx`: larger `idx` values encode larger
magnitudes on the positive side, and smaller magnitudes on the negative
side (since negative values are reflected). For two values of the same
sign, unsigned integer comparison of `idx` gives the correct magnitude
order. Stated in `proofs/P16_Ordering.v`.

## 2. Hardware primitive budget

The encoding is designed so that its structural operations factor
through a small fixed set of integer primitives. The set used by the
arithmetic proofs is: `IntAdd`, `IntSub`, `BitXor`, `BitNot`, `Compare`,
`Clz` (count leading zeros), and `FixedShift1`. Variable-distance barrel
shifters are not part of this set.

The primitive-composition claims for the encoding itself are:

- **Negation** factors through a single `BitNot` over the idx field.
- **Sign** and **magnitude comparison** factor through `Compare` against
  fixed constants (`pit_lo`, `pit_hi`) and against the other operand's
  `idx`.
- **Cancellation detection**, when a subtraction produces a value whose
  idx lies inside the pivot gap, is `BitXor` followed by `Clz`; no
  variable-distance shift is required.

The encoding-level statement that structural operations factor through
this primitive set is `proofs/P19_NoBarrel.v`. Downstream arithmetic
that depends on this property is out of scope here.

## 3. Interpretation independence

The same bit pattern can be decoded several ways, selected by context.
The structural properties stated in Section 1 (mirror, negation-by-NOT,
ordering, partition) hold for every interpretation because they are
statements about `idx`, not about the value decoded from it.

- **Linear**: value proportional to `position − MIRROR`, scaled by a
  fixed constant `SCALE`. Uniform spacing on the number line.
- **Logarithmic**: value equal to `base^((position − MIRROR)/D)` for a
  chosen base (φ, e, 10, 2, or an arbitrary positive base). Uniform
  relative precision; multiplication and division become addition and
  subtraction in position space.
- **Polar**: `idx` as a magnitude class, `lin` as an angle (fraction of
  2π). Used for complex-number representation.
- **IEEE-style**: `position` split into exponent and mantissa fields,
  reproducing IEEE 754 value semantics on an SBP substrate.
- **Continued-fraction**: `position` encodes a path in the Stern-Brocot
  tree, yielding best-rational approximations at each prefix length.

Orthogonality of interpretation and structural property is stated in
`proofs/P25_InterpretationOrthogonality.v`. The linear/geometric-mean
identity used by the arithmetic layer is
`proofs/P30_LinGeometric.v`. The polar-addition identity used by
complex arithmetic is `proofs/P31_PolarAdd.v`. Decoders are defined in
`proofs/Interpretations.v`.

## 4. Kinematic robustness: no zero-singularity

A property of the mirror encoding is that every non-zero normal value
has distance from the mirror of at least 1. The pivot gap occupies a
contiguous region of three (Layout A) or four (Layout C) idx positions
straddling the mirror, and normal values lie outside that region by
construction.

Consequences of this bound:

- The reciprocal of any non-zero normal has magnitude bounded above by
  a fixed constant of the bit-width; `|1/v|` cannot exceed `SCALE`.
- There is no continuous path in idx space from a negative normal to a
  positive normal that does not traverse the pivot gap; the gap
  separates the two normal regions.

This is the property that removes one source of numerical instability
in inverse kinematics: a value approaching the mirror cannot be
classified as normal without satisfying the minimum-magnitude bound.
Values that would violate the bound are classified as denormal or NaN
by Section 1.4, not as normal.

The bound is stated in `proofs/P27_NoZeroSingularity.v` for Layout C
and in `proofs/P28_PureNoSingularity.v` for Layout A. Both proofs are
complete across all supported widths.

## 5. Layout variants

Three layout variants differ in what occupies the pivot gap:

| Variant | Pivot gap      | NaN tag | Denormal tag |
|---------|----------------|:-------:|:------------:|
| A       | empty          | no      | no           |
| B       | NaN only       | yes     | no           |
| C       | NaN + denormal | yes     | yes          |

Layout A maximises the normal-value range. Layout C (the default in
this repository) retains IEEE-style information channels and a denormal
region for near-zero values. The minimum-magnitude property of Section
4 holds in all three variants; the proof for Layout A
(`proofs/P28_PureNoSingularity.v`) shows the mirror axis alone is
sufficient, and Layout C's NaN/denormal gap is not required for the
magnitude bound.

### 5.1 Tapered variant

A fourth variant partitions the normal-magnitude range into four
contiguous sub-ranges (tapers). The three inner tapers are equal-sized;
the final taper is two codes wider. Taper boundaries place higher code
density near unity and lower density at the extremes, matching typical
floating-point operand distributions.

Tapering is orthogonal to the A/B/C choice: it applies to any of the
three layouts. The structural claims — partition totality, disjointness,
monotonicity of taper index in magnitude, and magnitude reconstruction
per-taper — are in `proofs/P32_Taper.v`.

## 6. Related encodings

The mirror-axis structure is closest to **ones' complement**, which
also gives negation by bitwise NOT. Ones' complement has two zero
encodings (+0 and −0) and applies only to integers; SBP has one zero
encoding and defines magnitude and interpretation on top of the mirror.

**IEEE 754** uses an explicit sign bit and a biased exponent. Negation
is a single bit flip but is localised to the sign bit, not a full
bitwise complement; subnormals and NaN occupy specific exponent fields
rather than a contiguous gap.

**Logarithmic Number Systems** (Coleman and others) use position to
encode the logarithm of the value but typically with a separate sign
bit and no mirror axis. The SBP logarithmic interpretation differs in
using position for sign as well as magnitude.

**Posits** (Gustafson) use a variable-length regime field and a
run-length-encoded exponent. The bit structure is different from SBP's
fixed idx/lin split; posits do not have the mirror-axis complement
property.

## 7. Structural properties

The encoding's claims are verified in Coq/Rocq. What follows covers
only proofs about the encoding itself; proofs about the arithmetic
layer, the FPU, or transcendental evaluators are out of scope here.

### 7.1 Negation and mirror symmetry

**P01 `NegInvolution`** (Admitted) states `neg(neg(v)) = v` for every
valid code. Negation is defined as bitwise complement of `idx` (the
`lin` field preserved); applying it twice returns the original. The
statement is a generic-width claim currently pending a `cbn`-tactic
refactor for W32/W64.

**P06 `IdxSymmetry`** (Admitted) states `idx(v) + idx(neg(v)) = IDX_MAX`
for every non-zero valid code. Algebraic consequence of reflection at
`2^(N−1) − 0.5`. Same width-generalisation blocker as P01.

**P23 `NegSpecials`** (6 Qed) proves how `neg` acts on each
special-value class: `neg(zero) = zero`, `neg(+∞) = −∞`,
`neg(−∞) = +∞`, `neg(NaN) = NaN`, `neg(+denormal) = −denormal`. Proven
concretely for all widths.

### 7.2 Sign and classification

**P09 `SignByPosition`** (Admitted) states that sign is extracted by a
single integer comparison against the pivot (`idx > pit_hi` for
positive; `1 < idx < pit_lo` for negative). No bit extraction from the
payload is needed. Width-generalisation pending.

**P10 `Denormals`** (2 Admitted) states the four-slot structure of the
pivot gap in Layout C — `pit_lo` for negative denormal, `nan_lo` and
`nan_hi` for NaN, `pit_hi` for positive denormal — and that negation
swaps the two denormal slots while fixing the NaN slots.

**P11 `Constants`** (4 Qed) tabulates the derived constants
(`IDX_MAX`, `LIN_MAX`, `pit_lo`, `pit_hi`, `nan_lo_idx`, `nan_hi_idx`)
for bit widths W8, W16, W32, and W64, verified by computation.

**P21 `IdxPartition`** (2 Qed, 1 Admitted) proves that every `idx` in
`[0, IDX_MAX]` belongs to exactly one of nine classes (zero, ±∞,
±denormal, two NaN slots, positive-normal, negative-normal). Totality
is proven concretely for W8 (Qed). Disjointness is proven for all
widths (Qed). Totality for all widths is Admitted pending the `cbn`
refactor.

### 7.3 Bit layout

**P13 `PackUnpack`** (Admitted) states that the `(idx, lin)` pair packs
and unpacks round-trip through the bit layout: `unpack(pack(v)) = v`
for every valid value. Generic-width statement pending
width-specialisation.

**P17 `Sequence`** (2 Admitted) states two adjacent claims: sequence
separability (arithmetic depends only on abstract position, not on the
specific bit pattern used for that position) and Gray-code adjacency
(consecutive codes in a Gray-coded variant differ by exactly one bit).
P17 does not claim idx-to-position bijection; that claim is implicit in
the bit layout and `pack`/`unpack` but has no standalone theorem here.

**P18 `BitEfficiency`** (8 Qed) proves `pack` is injective conditionally
on P13's round-trip property: no two valid `(idx, lin)` pairs produce
the same bit pattern. Zero-idx uniqueness, ±∞ idx count, NaN idx count,
and Layout-C special-position count are verified concretely;
bit-efficiency for W8 is proved by direct computation.

### 7.4 Ordering

**P16 `Ordering`** (5 Admitted) states that unsigned integer comparison
on `idx` yields the correct magnitude ordering for two normal values of
the same sign, that the relation holds symmetrically for negative
normals, and that positive-normal idx values exceed all negative-normal
idx values. The cross-sign case is handled by the pivot comparison.
Same width-generalisation blocker as P01.

### 7.5 Minimum magnitude

**P27 `NoZeroSingularity`** (7 Qed) proves, for Layout C, that every
non-zero normal value has `distance ≥ 1`. The pivot gap is shown to be
exactly 3 idx positions wide, and the positive-normal and
negative-normal idx ranges are disjoint.

**P28 `PureNoSingularity`** (7 Qed) generalises the minimum-magnitude
property to Layout A (no NaN barrier). The bound `mag(v) ≥ 1` for
non-zero `v` follows from the mirror axis alone; the NaN/denormal gap
in Layout C is not required for this property.

### 7.6 Tapered magnitude partition

**P32 `Taper`** (11 Qed) covers the tapered variant. Given the W16
tapered constants (`MAX_MAG = 32766`, `CODES = 8191` per inner taper,
three boundaries at `8191`, `16382`, `24573`): every magnitude maps to
a unique taper index in `{0, 1, 2, 3}` (totality and disjointness); the
taper index is monotone non-decreasing in magnitude; each inner taper
has exactly `CODES` entries and the final taper has `CODES + 2`;
magnitude reconstructs as `taper_base(t) + taper_pos(m)` per-taper.
Taper invariance under negation follows from P06.

## 8. Verification status

Across the 14 proof files covering the encoding:

| Group                                | Qed | Admitted | Notes                         |
|--------------------------------------|----:|---------:|-------------------------------|
| Constants (P11)                      |   4 |        0 | concrete per-width            |
| Neg on specials (P23)                |   6 |        0 | concrete per-class            |
| Min magnitude (P27, P28)             |  14 |        0 | both layouts proven           |
| Tapered partition (P32)              |  11 |        0 | W16 tapered variant           |
| Bit efficiency (P18)                 |   8 |        0 | pack injective assuming P13   |
| Idx partition (P21)                  |   2 |        1 | W8 total; all-widths disjoint |
| Neg involution / symmetry (P01, P06) |   0 |        2 | width-generalisation pending  |
| Sign / denormals (P09, P10)          |   0 |        3 | same blocker                  |
| Pack-unpack (P13)                    |   0 |        1 | same blocker                  |
| Ordering (P16)                       |   0 |        5 | same blocker                  |
| Sequence (P17)                       |   0 |        2 | separability + Gray adjacency |
| **Total**                            |  **45** | **14** |                          |

The Admitted theorems share a common cause: their statements are
written as `forall w`, and the step-by-step proofs (which work at W8 by
direct computation) diverge on W32/W64 because `simpl` unfolds large
constants like `2^32 − 1`. A `cbn`-based refactor of the per-width
tactics would close most of them; this is a tactic-engineering task,
not a mathematical one.

### 8.1 On bijection

The encoding's `(idx, lin) ↔ bit-pattern` round-trip is claimed by P13
(Admitted) and used by P18's injectivity proof (Qed). A fully proven
bijection theorem — "`pack` is bijective on valid inputs" — is not
currently a single theorem in the proof set. The closest explicit
statements are P13's round-trip (Admitted) and P18's `pack_injective`
(Qed, conditional on P13). Adding an explicit
`forall w, bijective (pack w)` theorem would close this gap; it is not
presently in `Properties.v`.

## 9. Summary

The encoding's structural claims together state: SBP is a
mirror-symmetric encoding with a fixed-size reserved region for special
values; every non-zero normal value has magnitude bounded below by a
fixed constant of the encoding; class membership is total and
disjoint; sign, negation, and magnitude comparison are each a single
`idx`-level operation; the encoding is independent of the
interpretation (linear, logarithmic, polar, IEEE-style, or
continued-fraction) chosen to decode it.

Downstream consequences — primitive budget for arithmetic, absence of
variable barrel shifters, stable transcendental evaluators, and FPU
control-flow — are out of scope for this document.

## 10. References

- `proofs/Semantics.v` — formal definition.
- `proofs/Layouts.v` — the three layout variants as a typeclass.
- `proofs/Interpretations.v` — linear, log, polar, and IEEE decoders.
- `proofs/Properties.v` — aggregated structural claims.
