(** Ordering properties — integer compare gives correct float ordering *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

(** Positive values: idx monotonic with magnitude *)
Theorem idx_monotonic_positive_proof : forall w,
  idx_monotonic_positive w.
Admitted.

(** Negative values: idx anti-monotonic with magnitude *)
Theorem idx_monotonic_negative_proof : forall w,
  idx_monotonic_negative w.
Admitted.

(** Positive idx ordering iff magnitude ordering *)
Theorem positive_order_proof : forall w,
  positive_order_by_idx w.
Admitted.

(** Negative idx ordering iff magnitude ordering (reversed) *)
Theorem negative_order_proof : forall w,
  negative_order_by_idx w.
Admitted.

(** All positives above all negatives *)
Theorem positive_above_negative_proof : forall w,
  positive_above_negative w.
Admitted.
