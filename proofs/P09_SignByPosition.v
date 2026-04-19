(** Sign determined by position — the defining SBP property *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Properties.
Open Scope Z_scope.

Theorem sign_by_position_proof : forall w (v : Lsbp w),
  is_valid v -> sign_by_position_invariant v.
Admitted.
