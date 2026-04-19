(** * SBP Proof Tactics *)

From Stdlib Require Import ZArith.
From Stdlib Require Import Bool.
From Stdlib Require Import Lia.
Require Import Semantics.
Require Import Arithmetic.
Open Scope Z_scope.

(** Unfold all SBP definitions to bare Z conditionals *)
Ltac unfold_sbp :=
  unfold sign, lsbp_neg, lsbp_add, lsbp_mul,
    is_valid, is_nan, is_nan_b, is_zero, is_zero_b,
    is_inf, is_inf_b, is_pos_inf, is_pos_inf_b,
    is_neg_inf, is_neg_inf_b,
    is_positive, is_positive_b, is_negative, is_negative_b,
    is_normal, is_finite,
    is_pos_denormal, is_pos_denormal_b,
    is_neg_denormal, is_neg_denormal_b, is_denormal,
    denorm_pos_idx, denorm_neg_idx,
    distance, decode_pos, decode,
    lsbp_zero, lsbp_pos_inf, lsbp_neg_inf, lsbp_nan,
    lsbp_pos_denormal, lsbp_neg_denormal, lsbp_qnan,
    nan_lo_idx, nan_hi_idx,
    pit_hi, pit_lo, idx_max, lin_max, scale, half_bits
  in *.

(** One pass: resolve each Z.eqb/Z.ltb from left to right, max depth.
    Uses zstep for single resolution, zbool for bounded iteration. *)
Ltac zstep :=
  match goal with
  | |- context [Z.eqb ?a ?b] =>
      destruct (Z.eqb_spec a b);
      [subst; simpl; try reflexivity; try lia|simpl; try lia]
  | |- context [Z.ltb ?a ?b] =>
      destruct (Z.ltb_spec a b); simpl; try reflexivity; try lia
  | |- context [Z.leb ?a ?b] =>
      destruct (Z.leb_spec a b); simpl; try reflexivity; try lia
  | |- context [Bool.eqb ?a ?b] =>
      destruct a; destruct b; simpl; try reflexivity; try lia
  | |- context [andb ?a ?b] =>
      destruct a eqn:?; simpl; try reflexivity; try lia
  | |- context [orb ?a ?b] =>
      destruct a eqn:?; simpl; try reflexivity; try lia
  end.

(** Bounded resolution: try up to 12 steps then give up *)
Ltac zbool := try zstep; try zstep; try zstep; try zstep;
              try zstep; try zstep; try zstep; try zstep;
              try zstep; try zstep; try zstep; try zstep;
              try lia; try reflexivity.

(** Full tactic: unfold, destruct width, destruct record, resolve booleans *)
Ltac sbp_crush :=
  intros; unfold_sbp; simpl in *;
  match goal with
  | v : Lsbp _ |- _ => destruct v as [?i ?l]; simpl in *
  | _ => idtac
  end;
  zbool.
