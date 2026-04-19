(** * SBP Interpretations: Four ways to read a position

    SBP encoding maps values to positions. The position is an integer index.
    The INTERPRETATION defines what float each position represents.

    All four interpretations use the same SBP structure (Layout A/B/C).
    They differ only in the decode function: position -> real value.

    | Interpretation | decode(pos)              | mul cost | add cost  |
    |---------------|--------------------------|----------|-----------|
    | Linear        | sign * pos / scale       | O(n)     | O(1)      |
    | Log           | sign * base^(pos/D)      | O(1)     | O(log+exp)|
    | CF            | sign * cf_decode(pos)     | O(1)     | O(CF)     |
    | IEEE-style    | sign * 2^exp * (1 + frac) | O(n)     | O(1)      |
*)

From Stdlib Require Import ZArith.
From Stdlib Require Import QArith.
Require Import Semantics.
Open Scope Z_scope.

(** ** Linear Interpretation

    position maps linearly to value: val = sign * mag / SCALE
    This is what LSBP Semantics.v implements.

    Properties:
    - Uniform spacing: adjacent positions differ by 1/SCALE
    - Addition is position addition (free)
    - Multiplication requires full multiply (expensive)
    - Range limited by bit width: max value = idx_max / SCALE
*)

Record LinearInterp := mkLinearInterp {
  lin_scale : Z;   (* denominator: val = mag / lin_scale *)
}.

Definition linear_decode (interp : LinearInterp) (mag : Z) : Q :=
  Qmake mag (Z.to_pos (lin_scale interp)).

Definition linear_encode (interp : LinearInterp) (v : Q) : Z :=
  let num := Qnum v * lin_scale interp in
  let den := Z.pos (Qden v) in
  num / den.

(** Linear add: position add = value add *)
Definition linear_add (a b : Z) : Z := a + b.

(** Linear mul: requires full multiply + rescale *)
Definition linear_mul (interp : LinearInterp) (a b : Z) : Z :=
  (a * b) / lin_scale interp.

(** ** Logarithmic Interpretation

    position maps logarithmically: val = sign * base^((pos - ONE_POS) / D)
    where D = 2^lin_bits (subgrid resolution).

    Properties:
    - Exponential spacing: each step multiplies value by base^(1/D)
    - Multiplication is position addition (free): pos(a*b) = pos(a) + pos(b) - ONE_POS
    - Division is position subtraction (free): pos(a/b) = pos(a) - pos(b) + ONE_POS
    - Square root is position halving (free): pos(sqrt(a)) = (pos(a) + ONE_POS) / 2
    - Addition requires log-domain correction (expensive): log(a+b) = log(a) + log(1 + b/a)
*)

Record LogInterp := mkLogInterp {
  log_base_q48 : Z;  (* ln(base) in Q48 fixed point *)
  log_lin_bits : Z;   (* D = 2^lin_bits *)
  log_one_pos : Z;    (* position where value = 1.0 *)
}.

(** Log mul: just add positions, adjust for ONE_POS *)
Definition log_mul (interp : LogInterp) (pos_a pos_b : Z) : Z :=
  pos_a + pos_b - log_one_pos interp.

(** Log div: subtract positions, adjust for ONE_POS *)
Definition log_div (interp : LogInterp) (pos_a pos_b : Z) : Z :=
  pos_a - pos_b + log_one_pos interp.

(** Log sqrt: halve distance from ONE_POS *)
Definition log_sqrt (interp : LogInterp) (pos : Z) : Z :=
  (pos + log_one_pos interp) / 2.

(** Log add: requires base^(-gap/D) correction via CF or table *)
(** This is where the CF engine comes in — not formalized here,
    see the Rust parametric module for the implementation *)

(** Standard log bases *)
Definition log_phi : LogInterp := mkLogInterp 135132266347361024 0 0.
Definition log_e   : LogInterp := mkLogInterp 281474976710656 0 0.
Definition log_10  : LogInterp := mkLogInterp 648051823176498688 0 0.
Definition log_137 : LogInterp := mkLogInterp 1385334909009100800 0 0.
Definition log_128 : LogInterp := mkLogInterp 1365799573498470400 0 0.

(** ** CF (Continued Fraction) Interpretation

    position encodes a CF convergent p/q.
    Value = sign * p / q where p/q is the k-th convergent of some CF.

    Properties:
    - Best rational approximation for given denominator size
    - Multiplication: mediant or convergent composition
    - Angle representation: p/q fraction of 2*pi (polar)
    - Natural for Stern-Brocot tree traversal

    Used by CfPolar for angle encoding.
*)

Record CfInterp := mkCfInterp {
  cf_max_depth : Z;  (* maximum CF depth *)
}.

(** CF convergent from partial quotients.
    A position encodes a path in the Stern-Brocot tree:
    left = increment denominator, right = increment numerator.
    The position bits encode the path. *)

(** Mediant of two fractions: (a+c)/(b+d) *)
Definition mediant (p1 q1 p2 q2 : Z) : Z * Z :=
  (p1 + p2, q1 + q2).

(** ** IEEE-style Interpretation

    position splits into exponent and fraction: val = sign * 2^exp * (1 + frac/2^F)
    where F = fraction bits, exp = position >> F, frac = position mod 2^F.

    This is IEEE 754 semantics applied to SBP encoding.
    The key difference from actual IEEE: the ENCODING is SBP (mirror symmetry,
    sign by position), but the INTERPRETATION is IEEE (exponent + mantissa).

    Properties:
    - Variable spacing: spacing doubles with each exponent step
    - Addition requires alignment (shift to match exponents)
    - Multiplication: add exponents, multiply mantissas
    - Same precision trade-offs as IEEE, but better encoding properties
*)

Record IeeeInterp := mkIeeeInterp {
  ieee_exp_bits : Z;   (* number of exponent bits *)
  ieee_frac_bits : Z;  (* number of fraction/mantissa bits *)
  ieee_bias : Z;       (* exponent bias *)
}.

Definition ieee_exponent (interp : IeeeInterp) (pos : Z) : Z :=
  Z.shiftr pos (ieee_frac_bits interp) - ieee_bias interp.

Definition ieee_fraction (interp : IeeeInterp) (pos : Z) : Z :=
  Z.land pos (Z.shiftl 1 (ieee_frac_bits interp) - 1).

(** IEEE mul: add exponents, multiply fractions *)
Definition ieee_mul_exp (interp : IeeeInterp) (pos_a pos_b : Z) : Z :=
  ieee_exponent interp pos_a + ieee_exponent interp pos_b.

(** Standard IEEE formats as SBP interpretations *)
Definition ieee_half   : IeeeInterp := mkIeeeInterp 5 10 15.
Definition ieee_single : IeeeInterp := mkIeeeInterp 8 23 127.
Definition ieee_double : IeeeInterp := mkIeeeInterp 11 52 1023.

(** ** Interpretation Independence

    The four interpretations are ORTHOGONAL to SBP encoding.
    The same SBP bit pattern (Layout A/B/C) can be interpreted
    as linear, log, CF, or IEEE-style.

    SBP encoding properties (negation, ordering, CLZ) hold
    regardless of interpretation. The interpretation only affects
    what float value each position represents.

    | SBP property          | Depends on interpretation? |
    |-----------------------|---------------------------|
    | neg(neg(x)) = x       | No — pure encoding        |
    | sign by position       | No — pure encoding        |
    | integer comparison     | No — pure encoding        |
    | CLZ cancellation       | No — pure encoding        |
    | mul = pos add          | Yes — only for Log        |
    | add = pos add          | Yes — only for Linear     |
    | spacing uniformity     | Yes — Linear vs Log       |
    | dynamic range          | Yes — Log >> Linear       |
*)

(** The interpretation is a parameter, not a property of SBP *)
Inductive Interpretation :=
  | InterpLinear  : LinearInterp -> Interpretation
  | InterpLog     : LogInterp -> Interpretation
  | InterpCF      : CfInterp -> Interpretation
  | InterpIEEE    : IeeeInterp -> Interpretation.
