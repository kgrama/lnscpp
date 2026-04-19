(** Constants per width — verified by computation *)

From Stdlib Require Import ZArith.
Require Import Semantics.
Require Import Properties.
Open Scope Z_scope.

Theorem constants_8_proof : constants_8.
Proof. unfold constants_8; repeat split; reflexivity. Qed.

Theorem constants_16_proof : constants_16.
Proof. unfold constants_16; repeat split; reflexivity. Qed.

Theorem constants_32_proof : constants_32.
Proof. unfold constants_32; repeat split; reflexivity. Qed.

Theorem constants_64_proof : constants_64.
Proof. unfold constants_64; repeat split; reflexivity. Qed.
