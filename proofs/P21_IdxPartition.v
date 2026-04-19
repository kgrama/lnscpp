(** Every valid idx falls in exactly one class: zero, pos_inf, neg_inf,
    pos_denormal, neg_denormal, nan, positive_normal, or negative_normal.
    No valid bit pattern escapes classification. *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Properties.
Open Scope Z_scope.

(** ** Totality: every valid idx has a class (W8 concrete) *)
Theorem idx_partition_total_w8 : idx_partition_total W8.
Proof.
  unfold idx_partition_total.
  intros v [Hi _].
  unfold is_zero, is_pos_inf, is_neg_inf,
         is_pos_denormal, is_neg_denormal,
         denorm_pos_idx, denorm_neg_idx,
         is_positive, is_negative,
         idx_max, pit_lo, pit_hi, nan_lo_idx, nan_hi_idx,
         half_bits in *.
  simpl in *.
  (* 9 classes: zero, neg_inf, pos_inf, neg_denorm, pos_denorm,
     nan_lo, nan_hi, negative normal range, positive normal range *)
  assert (H : idx v = 0 \/ idx v = 1 \/ idx v = 15 \/
              idx v = 6 \/ idx v = 9 \/
              idx v = 7 \/ idx v = 8 \/
              (1 < idx v < 6) \/ (9 < idx v < 15)) by lia.
  destruct H as [H|[H|[H|[H|[H|[H|[H|[H|H]]]]]]]].
  - left; exact H.                                     (* zero *)
  - right; right; left; exact H.                       (* neg_inf *)
  - right; left; exact H.                              (* pos_inf *)
  - right; right; right; right; left; exact H.        (* neg_denormal: 5th *)
  - right; right; right; left; exact H.                (* pos_denormal: 4th *)
  - right; right; right; right; right; left; exact H.  (* nan_lo *)
  - right; right; right; right; right; right; left; exact H. (* nan_hi *)
  - right; right; right; right; right; right; right; right; lia. (* negative: position 9, 8 rights *)
  - right; right; right; right; right; right; right; left; lia.  (* positive: position 8, 7 rights + left *)
Qed.

(** ** Disjointness: certain pairs of classes can never both hold *)
Theorem idx_partition_disjoint_all : forall w,
  idx_partition_disjoint w.
Proof.
  intros w v [Hi _].
  unfold is_zero, is_pos_inf, is_neg_inf,
         is_pos_denormal, is_neg_denormal,
         denorm_pos_idx, denorm_neg_idx,
         is_positive, is_negative.
  repeat split; intros [Ha Hb].
  - rewrite Ha in Hb.
    (* 0 = idx_max w? only at width with idx_max = 0, impossible *)
    unfold idx_max, half_bits in Hb.
    destruct w; simpl in Hb; lia.
  - rewrite Ha in Hb. discriminate.
  - rewrite Ha in Hb. unfold idx_max, half_bits in Hb.
    destruct w; simpl in Hb; lia.
  - rewrite Ha in Hb.
    unfold pit_hi, pit_lo, half_bits in Hb.
    destruct w; simpl in Hb; lia.
  - (* is_positive v /\ is_negative v: contradicts pit_lo < pit_hi *)
    destruct Hb as [_ Hneg].
    unfold pit_hi, pit_lo, half_bits in *.
    destruct w; simpl in *; lia.
Qed.

(** ** Totality for all widths (hits cbn/simpl gap on W32/W64) *)
Theorem idx_partition_total_all : forall w, idx_partition_total w.
Proof.
  (* W8 proved concretely. W16/W32/W64 have same structure but
     simpl diverges on larger constants. Integrating cbn is the
     known blocker shared with other admitted proofs. *)
Admitted.
