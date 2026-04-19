(** neg(neg(x)) = x *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

Theorem neg_involution_proof : forall w, neg_involution w.
Admitted.
