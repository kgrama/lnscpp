(** Structural SBP operations (neg, sign, compare) are defined purely
    on Lsbp, independent of interpretation (Linear/Log/CF/IEEE).
    This justifies the "one FPU, mode-bit picks interpretation" design. *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

(** ** sign depends only on idx *)
Theorem sign_idx_only_proof : forall w, sign_idx_only w.
Proof.
  intros w a b Heq.
  unfold sign, is_pos_denormal_b, is_neg_denormal_b. rewrite Heq. reflexivity.
Qed.

(** ** Within positive half, idx order = distance order *)
Theorem compare_structural_proof : forall w, compare_structural w.
Proof.
  intros w a b Ha Hb Hpa Hpb.
  unfold distance, is_positive_b.
  unfold is_positive in Hpa, Hpb.
  assert (Ha' : pit_hi w <? idx a = true) by (apply Z.ltb_lt; lia).
  assert (Hb' : pit_hi w <? idx b = true) by (apply Z.ltb_lt; lia).
  rewrite Ha'. rewrite Hb'.
  repeat split; intros H; lia.
Qed.

(** ** neg's idx-level behavior is interpretation-free.
    For non-special values, neg(v).idx = idx_max - v.idx.
    For specials, neg maps inf<->inf, and interpretation is irrelevant. *)
Theorem neg_idx_interpretation_free_proof : forall w,
  neg_idx_interpretation_free w.
Proof.
  intros w v Hv Hnz Hnn Hnpd Hnnd.
  unfold lsbp_neg.
  destruct (Z.eqb_spec (idx v) 0) as [He | Hne].
  - unfold is_zero in Hnz. contradiction.
  - destruct (is_nan_b v) eqn:Enan.
    + exfalso. apply Hnn.
      unfold is_nan, is_nan_b in *.
      apply Bool.orb_true_iff in Enan.
      destruct Enan as [E | E];
        [left | right]; apply Z.eqb_eq; exact E.
    + destruct (Z.eqb_spec (idx v) (idx_max w)) as [Hm | Hnm].
      * right. left. split.
        -- unfold is_pos_inf. exact Hm.
        -- reflexivity.
      * destruct (Z.eqb_spec (idx v) 1) as [H1 | Hn1].
        -- right. right. split.
           ++ unfold is_neg_inf. exact H1.
           ++ reflexivity.
        -- destruct (is_pos_denormal_b v) eqn:Epd.
           ++ exfalso. apply Hnpd.
              unfold is_pos_denormal, is_pos_denormal_b in *.
              apply Z.eqb_eq. exact Epd.
           ++ destruct (is_neg_denormal_b v) eqn:End.
              ** exfalso. apply Hnnd.
                 unfold is_neg_denormal, is_neg_denormal_b in *.
                 apply Z.eqb_eq. exact End.
              ** left. simpl. reflexivity.
Qed.

(** ** Summary: structural ops are interpretation-orthogonal *)
Theorem sbp_structural_orthogonal : forall w,
  sign_idx_only w /\
  compare_structural w /\
  neg_idx_interpretation_free w.
Proof.
  intros w. repeat split.
  - intros a b H. apply sign_idx_only_proof; assumption.
  - apply compare_structural_proof; assumption.
  - apply compare_structural_proof; assumption.
  - apply compare_structural_proof; assumption.
  - intros v Hv Hnz Hnn Hnpd Hnnd.
    apply neg_idx_interpretation_free_proof; assumption.
Qed.
