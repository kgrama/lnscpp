(** Full specification of lsbp_neg on each special class.
    Complements P01 (involution) with class-by-class behavior. *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

(** ** neg(zero) = zero *)
Theorem neg_zero_proof : forall w, neg_zero_spec w.
Proof.
  intros w. unfold neg_zero_spec, lsbp_neg, lsbp_zero. reflexivity.
Qed.

(** ** neg(+inf) = -inf — works on W8/W16, uses cbv for larger widths *)
Theorem neg_pos_inf_proof : forall w, neg_pos_inf_spec w.
Proof.
  unfold neg_pos_inf_spec, lsbp_pos_inf, lsbp_neg_inf, lsbp_neg.
  intros w. destruct w; vm_compute; reflexivity.
Qed.

(** ** neg(-inf) = +inf *)
Theorem neg_neg_inf_proof : forall w, neg_neg_inf_spec w.
Proof.
  unfold neg_neg_inf_spec, lsbp_neg_inf, lsbp_pos_inf, lsbp_neg.
  intros w. destruct w; vm_compute; reflexivity.
Qed.

(** ** neg(+denormal) = -denormal *)
Theorem neg_pos_denormal_proof : forall w, neg_pos_denormal_spec w.
Proof.
  unfold neg_pos_denormal_spec, lsbp_pos_denormal, lsbp_neg_denormal, lsbp_neg.
  intros w l _. destruct w; vm_compute; reflexivity.
Qed.

(** ** neg(-denormal) = +denormal *)
Theorem neg_neg_denormal_proof : forall w, neg_neg_denormal_spec w.
Proof.
  unfold neg_neg_denormal_spec, lsbp_neg_denormal, lsbp_pos_denormal, lsbp_neg.
  intros w l _. destruct w; vm_compute; reflexivity.
Qed.

(** ** neg(NaN) is NaN *)
Theorem neg_nan_proof : forall w, neg_nan_stays_nan w.
Proof.
  intros w v Hnan. unfold lsbp_neg.
  destruct Hnan as [H | H].
  - (* idx v = nan_lo_idx w *)
    assert (Hnz : idx v <> 0).
    { rewrite H. unfold nan_lo_idx, half_bits. destruct w; simpl; lia. }
    destruct (Z.eqb_spec (idx v) 0) as [He | _]; [congruence|].
    unfold is_nan_b.
    assert (Heq : Z.eqb (idx v) (nan_lo_idx w) = true)
      by (apply Z.eqb_eq; exact H).
    rewrite Heq. simpl.
    left. exact H.
  - (* idx v = nan_hi_idx w *)
    assert (Hnz : idx v <> 0).
    { rewrite H. unfold nan_hi_idx, half_bits. destruct w; simpl; lia. }
    destruct (Z.eqb_spec (idx v) 0) as [He | _]; [congruence|].
    unfold is_nan_b.
    assert (Heq : Z.eqb (idx v) (nan_hi_idx w) = true)
      by (apply Z.eqb_eq; exact H).
    rewrite Heq. rewrite Bool.orb_true_r. simpl.
    right. exact H.
Qed.
