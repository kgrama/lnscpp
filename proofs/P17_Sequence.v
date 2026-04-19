(** Bit sequence separability — arithmetic is sequence-blind *)

From Stdlib Require Import ZArith Bool Lia.
Require Import Semantics.
Require Import Arithmetic.
Require Import Properties.
Open Scope Z_scope.

(** Arithmetic depends only on abstract position, not bit pattern *)
Theorem sequence_separable_proof : forall w,
  sequence_separable w.
Admitted.

(** Gray sequence adjacency: consecutive codes differ by 1 bit *)
Theorem gray_adjacency_proof : forall w,
  gray_adjacency w.
Admitted.
