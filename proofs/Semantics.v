(** * LSBP Semantics: Linear Sign-by-Position Encoding *)

(**
  LSBP uses pointer arithmetic: (idx, lin) encodes a value where:
  - idx = index into value space (sign by position relative to pivot)
  - lin = linear interpolation within the cell

  Pointer Arithmetic Model:
  - pos  = dist * LIN_MAX + lin   (combined representation)
  - mag  = pos - LIN_MAX          (scaled magnitude)
  - val  = sign * mag / SCALE     (decoded real value)

  Where dist = distance from pivot:
  - Positive: dist = idx - PIT_HI
  - Negative: dist = PIT_LO - idx

  Sign-by-Position:
  - idx > PIT_HI => positive
  - idx < PIT_LO => negative
  - idx in [PIT_LO, PIT_HI] => special (zero, NaN)

  Negation: neg(idx) = IDX_MAX - idx (bitwise complement)

  Width | idx bits | lin bits | IDX_MAX | LIN_MAX | SCALE
  ------|----------|----------|---------|---------|-------
    8   |    4     |    4     |    15   |    15   |   16
   16   |    8     |    8     |   255   |   255   |  256
   32   |   16     |   16     | 65535   | 65535   | 65536
   64   |   32     |   32     | 2^32-1  | 2^32-1  | 2^32
*)

From Stdlib Require Import ZArith.
From Stdlib Require Import QArith.
Open Scope Z_scope.

(** ** Bit Width Parameter *)

Inductive BitWidth := W8 | W16 | W32 | W64.

Definition half_bits (w : BitWidth) : Z :=
  match w with
  | W8  => 4
  | W16 => 8
  | W32 => 16
  | W64 => 32
  end.

(** ** Derived Constants *)

Definition idx_max (w : BitWidth) : Z := Z.shiftl 1 (half_bits w) - 1.
Definition lin_max (w : BitWidth) : Z := Z.shiftl 1 (half_bits w) - 1.
Definition scale (w : BitWidth) : Z := Z.shiftl 1 (half_bits w).

(** Pivot: center of idx space with gap for specials *)
Definition pit_lo (w : BitWidth) : Z := Z.shiftl 1 (half_bits w - 1) - 2.
Definition pit_hi (w : BitWidth) : Z := Z.shiftl 1 (half_bits w - 1) + 1.

(** Special indices within pivot gap [PIT_LO..PIT_HI]:
    PIT_LO = denorm_neg (outer)
    NAN_LO = nan (inner)
    NAN_HI = nan (inner)
    PIT_HI = denorm_pos (outer)
*)
Definition nan_lo_idx (w : BitWidth) : Z := Z.shiftl 1 (half_bits w - 1) - 1.
Definition nan_hi_idx (w : BitWidth) : Z := Z.shiftl 1 (half_bits w - 1).
Definition denorm_pos_idx (w : BitWidth) : Z := pit_hi w.
Definition denorm_neg_idx (w : BitWidth) : Z := pit_lo w.

(** ** LSBP Value Representation *)

Record Lsbp (w : BitWidth) := mkLsbp {
  idx : Z;
  lin : Z;
}.

Arguments idx {w}.
Arguments lin {w}.
Arguments mkLsbp {w}.

(** ** Classification Predicates *)

Definition is_valid {w} (v : Lsbp w) : Prop :=
  0 <= idx v <= idx_max w /\ 0 <= lin v <= lin_max w.

Definition is_zero {w} (v : Lsbp w) : Prop := idx v = 0.
Definition is_pos_inf {w} (v : Lsbp w) : Prop := idx v = idx_max w.
Definition is_neg_inf {w} (v : Lsbp w) : Prop := idx v = 1.
Definition is_inf {w} (v : Lsbp w) : Prop := is_pos_inf v \/ is_neg_inf v.

Definition is_nan {w} (v : Lsbp w) : Prop :=
  idx v = nan_lo_idx w \/ idx v = nan_hi_idx w.

(** Denormals: tiny values using outer edges of pivot gap *)
Definition is_pos_denormal {w} (v : Lsbp w) : Prop := idx v = denorm_pos_idx w.
Definition is_neg_denormal {w} (v : Lsbp w) : Prop := idx v = denorm_neg_idx w.
Definition is_denormal {w} (v : Lsbp w) : Prop := is_pos_denormal v \/ is_neg_denormal v.

Definition is_positive {w} (v : Lsbp w) : Prop := idx v > pit_hi w.
Definition is_negative {w} (v : Lsbp w) : Prop := idx v > 1 /\ idx v < pit_lo w.
Definition is_normal {w} (v : Lsbp w) : Prop := is_positive v \/ is_negative v.
Definition is_finite {w} (v : Lsbp w) : Prop := ~is_inf v /\ ~is_nan v.

(** ** Boolean Predicates *)

Definition is_zero_b {w} (v : Lsbp w) : bool := Z.eqb (idx v) 0.
Definition is_pos_inf_b {w} (v : Lsbp w) : bool := Z.eqb (idx v) (idx_max w).
Definition is_neg_inf_b {w} (v : Lsbp w) : bool := Z.eqb (idx v) 1.
Definition is_inf_b {w} (v : Lsbp w) : bool := is_pos_inf_b v || is_neg_inf_b v.

Definition is_nan_b {w} (v : Lsbp w) : bool :=
  Z.eqb (idx v) (nan_lo_idx w) || Z.eqb (idx v) (nan_hi_idx w).

Definition is_pos_denormal_b {w} (v : Lsbp w) : bool := Z.eqb (idx v) (denorm_pos_idx w).
Definition is_neg_denormal_b {w} (v : Lsbp w) : bool := Z.eqb (idx v) (denorm_neg_idx w).
Definition is_denormal_b {w} (v : Lsbp w) : bool := is_pos_denormal_b v || is_neg_denormal_b v.

Definition is_positive_b {w} (v : Lsbp w) : bool := Z.ltb (pit_hi w) (idx v).
Definition is_negative_b {w} (v : Lsbp w) : bool :=
  Z.ltb 1 (idx v) && Z.ltb (idx v) (pit_lo w).

(** ** Sign Extraction *)

Definition sign {w} (v : Lsbp w) : Z :=
  if Z.eqb (idx v) 0 then 0
  else if Z.ltb (pit_hi w) (idx v) then 1              (* positive normal *)
  else if is_pos_denormal_b v then 1                   (* positive denormal *)
  else if Z.ltb 1 (idx v) && Z.ltb (idx v) (pit_lo w) then (-1)  (* negative normal *)
  else if is_neg_denormal_b v then (-1)                (* negative denormal *)
  else 0.                                              (* zero/NaN *)

(** ** Distance from Pivot *)

Definition distance {w} (v : Lsbp w) : Z :=
  if is_positive_b v then idx v - pit_hi w
  else if is_negative_b v then pit_lo w - idx v
  else 0.

(** ** Special Values *)

Definition lsbp_zero (w : BitWidth) : Lsbp w := mkLsbp 0 0.
Definition lsbp_pos_inf (w : BitWidth) : Lsbp w := mkLsbp (idx_max w) 0.
Definition lsbp_neg_inf (w : BitWidth) : Lsbp w := mkLsbp 1 0.
Definition lsbp_nan (w : BitWidth) : Lsbp w := mkLsbp (nan_lo_idx w) 0.
Definition lsbp_pos_denormal {w} (l : Z) : Lsbp w := mkLsbp (denorm_pos_idx w) l.
Definition lsbp_neg_denormal {w} (l : Z) : Lsbp w := mkLsbp (denorm_neg_idx w) l.

Definition lsbp_qnan {w} (payload : Z) : Lsbp w :=
  mkLsbp (nan_lo_idx w) (Z.lor payload (Z.shiftl 1 (half_bits w - 1))).

(** ** Decode *)

Definition decode_pos {w} (v : Lsbp w) : Z :=
  distance v * lin_max w + lin v.

(** Denormal decoding: val = sign * lin / SCALE *)
Definition decode_denormal {w} (v : Lsbp w) : Q :=
  let s := scale w in
  Qmult (inject_Z (sign v)) (Qmake (lin v) (Z.to_pos s)).

Definition decode {w} (v : Lsbp w) : Q :=
  if is_zero_b v then 0
  else if is_nan_b v then 0
  else if is_denormal_b v then decode_denormal v
  else
    let pos := decode_pos v - lin_max w in
    let s := scale w in
    Qmult (inject_Z (sign v)) (Qmake pos (Z.to_pos s)).

(** ** Encode *)

Definition encode_from_pos {w} (pos : Z) (positive : bool) : Lsbp w :=
  let lm := lin_max w in
  let adjusted := pos + lm in
  let d := adjusted / lm in
  let l := adjusted mod lm in
  if positive
  then mkLsbp (pit_hi w + d) l
  else mkLsbp (pit_lo w - d) l.

(** ** Pivot Gap *)

Definition in_pivot_gap {w} (i : Z) : Prop :=
  pit_lo w <= i <= pit_hi w.

(** ** Well-formedness *)

Definition wf_lsbp {w} (v : Lsbp w) : Prop :=
  is_valid v /\
  (is_zero v -> lin v = 0) /\
  (is_inf v -> lin v = 0).

(** ** Width Scaling *)

(**
  Scale factor between widths: 2^(half_bits(w2) - half_bits(w1))
  Scaling preserves value without re-encoding through decode/encode.
*)

Definition scale_factor (w1 w2 : BitWidth) : Z :=
  Z.shiftl 1 (half_bits w2 - half_bits w1).

(** ** Width Scaling (Pointer Arithmetic)

    Key insight: pos = dist * LIN_MAX + lin
    When scaling: pos2 = pos1 * scale_factor
    Since LIN_MAX scales proportionally, dist stays same, only lin scales.

    Scale up:   dist unchanged, lin *= sf
    Scale down: dist unchanged, lin /= sf
*)

(** Scale up: W8 -> W16 -> W32 -> W64 *)
Definition scale_up (w1 w2 : BitWidth) (v : Lsbp w1) : Lsbp w2 :=
  let sf := scale_factor w1 w2 in
  if is_zero_b v then lsbp_zero w2
  else if is_nan_b v then lsbp_nan w2
  else if is_pos_inf_b v then lsbp_pos_inf w2
  else if is_neg_inf_b v then lsbp_neg_inf w2
  else if is_pos_denormal_b v then lsbp_pos_denormal (lin v * sf)
  else if is_neg_denormal_b v then lsbp_neg_denormal (lin v * sf)
  else if is_positive_b v then
    let dist := idx v - pit_hi w1 in
    mkLsbp (pit_hi w2 + dist) (lin v * sf)
  else (* negative *)
    let dist := pit_lo w1 - idx v in
    mkLsbp (pit_lo w2 - dist) (lin v * sf).

(** Scale down: W64 -> W32 -> W16 -> W8 *)
Definition scale_down (w1 w2 : BitWidth) (v : Lsbp w1) : Lsbp w2 :=
  let sf := scale_factor w2 w1 in  (* inverse: w2 is smaller *)
  if is_zero_b v then lsbp_zero w2
  else if is_nan_b v then lsbp_nan w2
  else if is_pos_inf_b v then lsbp_pos_inf w2
  else if is_neg_inf_b v then lsbp_neg_inf w2
  else if is_pos_denormal_b v then lsbp_pos_denormal (lin v / sf)
  else if is_neg_denormal_b v then lsbp_neg_denormal (lin v / sf)
  else if is_positive_b v then
    let dist := idx v - pit_hi w1 in
    mkLsbp (pit_hi w2 + dist) (lin v / sf)
  else (* negative *)
    let dist := pit_lo w1 - idx v in
    mkLsbp (pit_lo w2 - dist) (lin v / sf).

(** Direct bit reinterpretation for 2x scaling *)
Definition pack {w} (v : Lsbp w) : Z :=
  idx v * Z.shiftl 1 (half_bits w) + lin v.

Definition unpack (w : BitWidth) (bits : Z) : Lsbp w :=
  let hb := half_bits w in
  let mask := Z.shiftl 1 hb - 1 in
  mkLsbp (Z.shiftr bits hb) (Z.land bits mask).
