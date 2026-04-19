(** Hardware primitive budget: SBP structural operations never need
    a variable-distance shifter (barrel). Each structural op factors
    through HwPrim = {IntAdd, IntSub, BitXor, BitNot, Compare, Clz,
    FixedShift1}. *)

From Stdlib Require Import ZArith List.
Require Import Semantics.
Require Import Properties.
Import ListNotations.
Open Scope Z_scope.

(** ** neg factors through BitNot alone *)
Theorem neg_no_barrel :
  op_uses_only_hwprim neg_circuit.
Proof.
  unfold op_uses_only_hwprim, neg_circuit.
  intros p Hin. simpl in Hin.
  destruct Hin as [He | []]. subst. right; right; right; left. reflexivity.
Qed.

(** ** sign factors through Compare alone *)
Theorem sign_no_barrel :
  op_uses_only_hwprim sign_circuit.
Proof.
  unfold op_uses_only_hwprim, sign_circuit.
  intros p Hin. simpl in Hin.
  destruct Hin as [He | []]. subst. right; right; right; right; left. reflexivity.
Qed.

(** ** compare factors through Compare alone *)
Theorem compare_no_barrel :
  op_uses_only_hwprim compare_circuit.
Proof.
  unfold op_uses_only_hwprim, compare_circuit.
  intros p Hin. simpl in Hin.
  destruct Hin as [He | []]. subst. right; right; right; right; left. reflexivity.
Qed.

(** ** classify factors through Compare alone *)
Theorem classify_no_barrel :
  op_uses_only_hwprim classify_circuit.
Proof.
  unfold op_uses_only_hwprim, classify_circuit.
  intros p Hin. simpl in Hin.
  destruct Hin as [He | []]. subst. right; right; right; right; left. reflexivity.
Qed.

(** ** cancellation detection uses BitXor + Clz, no variable shift *)
Theorem cancel_detect_no_barrel :
  op_uses_only_hwprim cancel_detect_circuit.
Proof.
  unfold op_uses_only_hwprim, cancel_detect_circuit.
  intros p Hin. simpl in Hin.
  destruct Hin as [He | [He | []]]; subst.
  - right; right; left. reflexivity.
  - right; right; right; right; right; left. reflexivity.
Qed.

(** ** Summary: no structural op requires a barrel shifter *)
Theorem sbp_no_barrel_summary :
  op_uses_only_hwprim neg_circuit /\
  op_uses_only_hwprim sign_circuit /\
  op_uses_only_hwprim compare_circuit /\
  op_uses_only_hwprim classify_circuit /\
  op_uses_only_hwprim cancel_detect_circuit.
Proof.
  repeat split;
    [ apply neg_no_barrel
    | apply sign_no_barrel
    | apply compare_no_barrel
    | apply classify_no_barrel
    | apply cancel_detect_no_barrel ].
Qed.
