(** Denormal sign invariant and negation swaps *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

Theorem denormal_sign_proof : forall w (v : Lsbp w),
  is_valid v -> denormal_sign_invariant v.
Admitted.

Theorem denormal_neg_swaps_proof : forall w, denormal_neg_swaps w.
Admitted.
