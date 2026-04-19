(** pack/unpack round-trip *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Properties.
Open Scope Z_scope.

Theorem pack_unpack_proof : forall w, pack_unpack_identity w.
Admitted.
