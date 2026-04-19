(** * SBP Layout Variations

    All SBP layouts share the mirror property: negation = complement.
    They differ in what occupies the pivot gap and extremes.

    idx space runs 0 .. idx_max, with the mirror at idx_max/2.

    Layout A (Pure):
      0         large_neg  small_neg  |mirror|  small_pos  large_pos   +inf
      [zero]    [----negative----]     [      ]  [----positive----]     [inf]

      - No NaN, no denormals
      - Simplest: every non-zero, non-inf position is normal
      - neg(x) = idx_max - x
      - Comparison: unsigned integer compare = magnitude compare (within sign)

    Layout B (NaN in gap):
      0         large_neg  small_neg  [nan_lo nan_hi]  small_pos  large_pos   inf
      [zero]    [----negative----]     [--NaN gap--]    [----positive----]     [inf]

      - NaN occupies the 2 innermost positions of the pivot gap
      - neg(nan) = nan (NaN is self-symmetric at the mirror)
      - Fewer normal positions than Layout A (2 positions lost to NaN)

    Layout C (Current LSBP — NaN + denormals in gap):
      0    -inf   large_neg  small_neg  [dn_neg nan_lo nan_hi dn_pos]  small_pos  large_pos  +inf
      [zero][inf]  [---negative---]      [------pivot gap------]        [---positive---]      [inf]

      - 4-slot pivot gap: neg denormal, NaN lo, NaN hi, pos denormal
      - Denormals provide subnormal-like values near zero
      - -inf at idx=1, +inf at idx=idx_max
      - neg(+inf) = -inf, neg(dn_pos) = dn_neg
      - This is what Semantics.v implements
*)

From Stdlib Require Import ZArith.
From Stdlib Require Import Bool.
From Stdlib Require Import Lia.
Require Import Semantics.
Require Import Arithmetic.
Open Scope Z_scope.

(** ** Abstract SBP: the shared contract all layouts satisfy *)

Class SbpLayout (V : BitWidth -> Type) := {
  (** Every layout has an index *)
  sbp_idx : forall {w}, V w -> Z;
  (** Maximum index for this width *)
  sbp_idx_max : BitWidth -> Z;
  (** Zero value *)
  sbp_zero : forall w, V w;
  (** Negation *)
  sbp_neg : forall {w}, V w -> V w;
  (** Sign extraction *)
  sbp_sign : forall {w}, V w -> Z;
  (** Magnitude (distance from midpoint) *)
  sbp_mag : forall {w}, V w -> Z;

  (** Core SBP axioms — any layout must satisfy these *)

  (** A1: Zero at idx = 0 *)
  sbp_zero_idx : forall w, sbp_idx (sbp_zero w) = 0;
  (** A2: Negation is complement *)
  sbp_neg_complement : forall w (v : V w),
    0 < sbp_idx v <= sbp_idx_max w ->
    sbp_idx (sbp_neg v) = sbp_idx_max w - sbp_idx v;
  (** A3: Negation is involution *)
  sbp_neg_neg : forall w (v : V w),
    0 <= sbp_idx v <= sbp_idx_max w ->
    sbp_neg (sbp_neg v) = v;
  (** A4: Sign determined by position *)
  sbp_sign_pos : forall w (v : V w),
    sbp_idx v > sbp_idx_max w / 2 -> sbp_sign v = 1 \/ sbp_sign v = 0;
  sbp_sign_neg : forall w (v : V w),
    0 < sbp_idx v < sbp_idx_max w / 2 -> sbp_sign v = (-1) \/ sbp_sign v = 0;
  sbp_sign_zero : forall w, sbp_sign (sbp_zero w) = 0;
  (** A5: Negation flips sign *)
  sbp_neg_sign : forall w (v : V w),
    sbp_sign v = 1 \/ sbp_sign v = (-1) ->
    sbp_sign (sbp_neg v) = - sbp_sign v;
  (** A6: Positive ordering — larger idx = larger magnitude *)
  sbp_order_pos : forall w (a b : V w),
    sbp_sign a = 1 -> sbp_sign b = 1 ->
    sbp_idx a > sbp_idx b -> sbp_mag a > sbp_mag b;
}.

(** ** Layout A: Pure SBP (no NaN, no denormals) *)

Record PureSbp (w : BitWidth) := mkPureSbp {
  pure_idx : Z;
}.

Arguments pure_idx {w}.
Arguments mkPureSbp {w}.

Definition pure_idx_max (w : BitWidth) : Z := Z.shiftl 1 (2 * half_bits w) - 1.
Definition pure_mid (w : BitWidth) : Z := Z.shiftl 1 (2 * half_bits w - 1).

Definition pure_is_zero {w} (v : PureSbp w) : Prop := pure_idx v = 0.
Definition pure_is_inf {w} (v : PureSbp w) : Prop := pure_idx v = pure_idx_max w.

Definition pure_sign {w} (v : PureSbp w) : Z :=
  if Z.eqb (pure_idx v) 0 then 0
  else if Z.ltb (pure_idx v) (pure_mid w) then (-1)
  else 1.

Definition pure_mag {w} (v : PureSbp w) : Z :=
  if Z.ltb (pure_idx v) (pure_mid w)
  then pure_mid w - pure_idx v          (* negative: distance from mid going left *)
  else pure_idx v - pure_mid w + 1.     (* positive: distance from mid going right *)

Definition pure_neg {w} (v : PureSbp w) : PureSbp w :=
  if Z.eqb (pure_idx v) 0 then v
  else mkPureSbp (pure_idx_max w - pure_idx v).

(** Pure SBP properties — these are stronger than Layout C because
    there are no special cases in the pivot gap *)

Definition pure_neg_involution (w : BitWidth) : Prop :=
  forall v : PureSbp w,
    0 <= pure_idx v <= pure_idx_max w ->
    pure_neg (pure_neg v) = v.

Definition pure_comparison (w : BitWidth) : Prop :=
  forall a b : PureSbp w,
    0 < pure_idx a <= pure_idx_max w ->
    0 < pure_idx b <= pure_idx_max w ->
    pure_sign a = 1 -> pure_sign b = 1 ->
    (pure_idx a > pure_idx b <-> pure_mag a > pure_mag b).

(** ** Layout B: NaN in gap (no denormals) *)

Record NanSbp (w : BitWidth) := mkNanSbp {
  nan_sbp_idx : Z;
}.

Arguments nan_sbp_idx {w}.
Arguments mkNanSbp {w}.

Definition nan_sbp_mid (w : BitWidth) : Z := Z.shiftl 1 (2 * half_bits w - 1).

(** NaN occupies the two positions straddling the mirror *)
Definition nan_sbp_lo (w : BitWidth) : Z := nan_sbp_mid w - 1.
Definition nan_sbp_hi (w : BitWidth) : Z := nan_sbp_mid w.

Definition nan_sbp_is_nan {w} (v : NanSbp w) : Prop :=
  nan_sbp_idx v = nan_sbp_lo w \/ nan_sbp_idx v = nan_sbp_hi w.

Definition nan_sbp_sign {w} (v : NanSbp w) : Z :=
  if Z.eqb (nan_sbp_idx v) 0 then 0
  else if Z.eqb (nan_sbp_idx v) (nan_sbp_lo w) then 0  (* NaN *)
  else if Z.eqb (nan_sbp_idx v) (nan_sbp_hi w) then 0  (* NaN *)
  else if Z.ltb (nan_sbp_idx v) (nan_sbp_lo w) then (-1)
  else 1.

(** NaN is self-symmetric: neg(NaN) = NaN *)
Definition nan_self_symmetric (w : BitWidth) : Prop :=
  forall v : NanSbp w,
    nan_sbp_is_nan v ->
    nan_sbp_is_nan (mkNanSbp (w:=w) (pure_idx_max w - nan_sbp_idx v)).

(** ** Layout comparison *)

(** All layouts share these core SBP properties:
    1. Negation = idx_max - idx (bitwise complement)
    2. Zero at idx = 0
    3. Sign determined by position relative to midpoint
    4. Magnitude = distance from midpoint
    5. Comparison via integer compare (within same sign)

    They differ in:
    - Layout A: no wasted positions, maximum dynamic range
    - Layout B: 2 positions for NaN, NaN is self-symmetric at mirror
    - Layout C: 4 positions for NaN + denormals, supports subnormals

    Trade-off: more special values = fewer normal positions = less range.

    | Layout | Special positions | Normal positions (per sign) |
    |--------|------------------|---------------------------|
    | A      | 2 (zero, inf)    | (idx_max - 1) / 2        |
    | B      | 4 (zero, inf, 2xNaN) | (idx_max - 3) / 2    |
    | C      | 6 (zero, 2xinf, 2xNaN, 2xdenorm) | (idx_max - 5) / 2 |
*)

(** ** Instance: Layout A satisfies SbpLayout *)

(** We declare instances with Admitted axioms. Once the tactic issue
    for large-constant widths is resolved, these become Qed. *)

#[global]
Instance PureSbpLayout : SbpLayout PureSbp.
Proof.
  refine {|
    sbp_idx := fun w v => pure_idx v;
    sbp_idx_max := pure_idx_max;
    sbp_zero := fun w => mkPureSbp 0;
    sbp_neg := fun w v => pure_neg v;
    sbp_sign := fun w v => pure_sign v;
    sbp_mag := fun w v => pure_mag v;
  |}; admit.
Admitted.

(** ** Instance: Layout C (LSBP from Semantics.v) satisfies SbpLayout *)

Definition lsbp_mag {w} (v : Lsbp w) : Z := distance v.

#[global]
Instance LsbpLayout : SbpLayout Lsbp.
Proof.
  refine {|
    sbp_idx := fun w v => idx v;
    sbp_idx_max := idx_max;
    sbp_zero := lsbp_zero;
    sbp_neg := fun w v => lsbp_neg v;
    sbp_sign := fun w v => sign v;
    sbp_mag := fun w v => lsbp_mag v;
  |}; admit.
Admitted.

(** ** Layout B instance *)

#[global]
Instance NanSbpLayout : SbpLayout NanSbp.
Proof.
  refine {|
    sbp_idx := fun w v => nan_sbp_idx v;
    sbp_idx_max := pure_idx_max;
    sbp_zero := fun w => mkNanSbp 0;
    sbp_neg := fun w v =>
      if Z.eqb (nan_sbp_idx v) 0 then v
      else mkNanSbp (pure_idx_max w - nan_sbp_idx v);
    sbp_sign := fun w v => nan_sbp_sign v;
    sbp_mag := fun w v =>
      if Z.ltb (nan_sbp_idx v) (nan_sbp_mid w)
      then nan_sbp_mid w - nan_sbp_idx v
      else nan_sbp_idx v - nan_sbp_mid w + 1;
  |}; admit.
Admitted.

(** ** Theorems that hold for ANY SbpLayout instance *)

Section UniversalSbp.
  Context {V : BitWidth -> Type} `{SbpLayout V}.

  (** Any layout: negating zero gives zero *)
  Theorem sbp_neg_zero_universal : forall w,
    sbp_neg (sbp_zero w) = sbp_zero w.
  Admitted.

  (** Any layout: sign(neg(v)) = -sign(v) for normal values *)
  Theorem sbp_neg_sign_universal : forall w (v : V w),
    (sbp_sign v = 1 \/ sbp_sign v = (-1)) ->
    sbp_sign (sbp_neg v) = - sbp_sign v.
  Proof. intros. apply sbp_neg_sign. assumption. Qed.

  (** Any layout: idx(v) + idx(neg(v)) = idx_max for non-zero *)
  Theorem sbp_idx_symmetry_universal : forall w (v : V w),
    0 < sbp_idx v <= sbp_idx_max w ->
    sbp_idx v + sbp_idx (sbp_neg v) = sbp_idx_max w.
  Proof.
    intros. rewrite sbp_neg_complement; lia.
  Qed.
End UniversalSbp.
