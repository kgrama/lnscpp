(** idx(v) + idx(neg(v)) = idx_max *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

Theorem idx_symmetry_proof : forall w, idx_symmetry w.
Admitted.
