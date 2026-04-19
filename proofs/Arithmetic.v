(** * LSBP Arithmetic: Generic Operations *)

From Stdlib Require Import ZArith.
From Stdlib Require Import QArith.
From Stdlib Require Import Bool.
Require Import Semantics.
Open Scope Z_scope.

(** ** Negation (Complement)

    neg(idx) = IDX_MAX - idx (equivalent to bitwise ~idx)
    This flips the position relative to the center of idx space.
    lin is preserved (interpolation doesn't change sign).
*)

Definition lsbp_neg {w} (v : Lsbp w) : Lsbp w :=
  let i := idx v in
  let im := idx_max w in
  if Z.eqb i 0 then v                    (* zero stays zero *)
  else if is_nan_b v then v              (* NaN stays NaN *)
  else if Z.eqb i im then lsbp_neg_inf w (* +inf -> -inf *)
  else if Z.eqb i 1 then lsbp_pos_inf w  (* -inf -> +inf *)
  else if is_pos_denormal_b v then lsbp_neg_denormal (lin v)  (* +denorm -> -denorm *)
  else if is_neg_denormal_b v then lsbp_pos_denormal (lin v)  (* -denorm -> +denorm *)
  else mkLsbp (im - i) (lin v).          (* complement idx *)

Definition neg_normal {w} (v : Lsbp w) : Lsbp w :=
  mkLsbp (idx_max w - idx v) (lin v).

(** ** Addition Helpers

    Pointer arithmetic model:
    - pos = dist * LIN_MAX + lin  (combined representation)
    - mag = pos - LIN_MAX         (actual scaled magnitude)
    - val = sign * mag / SCALE    (decoded value)

    Addition: compute magnitudes, add/subtract, re-encode.
*)

Definition add_same_sign {w} (ap bp : Z) (positive : bool) : Lsbp w :=
  let lm := lin_max w in
  let sum := ap + bp + lm in
  let max_pos := pit_lo w * lm in
  if Z.leb max_pos sum then
    if positive then lsbp_pos_inf w else lsbp_neg_inf w
  else
    let rd := sum / lm in
    let rl := sum mod lm in
    if positive
    then mkLsbp (pit_hi w + rd) rl
    else mkLsbp (pit_lo w - rd) rl.

Definition add_diff_sign {w} (ap bp : Z) (a_pos : bool) : Lsbp w :=
  let lm := lin_max w in
  if Z.eqb ap bp then lsbp_zero w
  else
    let '(diff, res_pos) := if Z.ltb bp ap then (ap - bp, a_pos) else (bp - ap, negb a_pos) in
    let result := diff + lm in
    let rd := result / lm in
    let rl := result mod lm in
    if Z.eqb rd 0 then lsbp_zero w
    else if res_pos
    then mkLsbp (pit_hi w + rd) rl
    else mkLsbp (pit_lo w - rd) rl.

(** ** Addition *)

Definition lsbp_add {w} (a b : Lsbp w) : Lsbp w :=
  if is_nan_b a then lsbp_qnan (lin a)
  else if is_nan_b b then lsbp_qnan (lin b)
  else if is_zero_b a then b
  else if is_zero_b b then a
  else if is_inf_b a && is_inf_b b then
    if Bool.eqb (is_pos_inf_b a) (is_pos_inf_b b) then a
    else lsbp_nan w
  else if is_inf_b a then a
  else if is_inf_b b then b
  else
    let lm := lin_max w in
    let a_pos := is_positive_b a in
    let b_pos := is_positive_b b in
    let ad := distance a in
    let bd := distance b in
    let ap := ad * lm + lin a - lm in
    let bp := bd * lm + lin b - lm in
    if Bool.eqb a_pos b_pos
    then add_same_sign ap bp a_pos
    else add_diff_sign ap bp a_pos.

(** ** Subtraction *)

Definition lsbp_sub {w} (a b : Lsbp w) : Lsbp w :=
  lsbp_add a (lsbp_neg b).

(** ** Multiplication *)

Definition lsbp_mul {w} (a b : Lsbp w) : Lsbp w :=
  if is_nan_b a then lsbp_qnan (lin a)
  else if is_nan_b b then lsbp_qnan (lin b)
  else if (is_zero_b a && is_inf_b b) || (is_zero_b b && is_inf_b a) then lsbp_nan w
  else if is_zero_b a || is_zero_b b then lsbp_zero w
  else
    let lm := lin_max w in
    let sc := scale w in
    let a_pos := is_positive_b a in
    let b_pos := is_positive_b b in
    let res_pos := Bool.eqb a_pos b_pos in
    if is_inf_b a || is_inf_b b then
      if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
    else
      let ad := distance a in
      let bd := distance b in
      let ap := (ad * lm + lin a) - lm in
      let bp := (bd * lm + lin b) - lm in
      let full_prod := ap * bp in
      let prod := full_prod / sc + lm in
      let max_pos := pit_lo w * lm in
      if Z.leb max_pos prod then
        if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
      else
        let rd := prod / lm in
        let rl := prod mod lm in
        if Z.eqb rd 0 then lsbp_zero w
        else if res_pos
        then mkLsbp (pit_hi w + rd) rl
        else mkLsbp (pit_lo w - rd) rl.

(** ** Division *)

Definition lsbp_div {w} (a b : Lsbp w) : Lsbp w :=
  if is_nan_b a then lsbp_qnan (lin a)
  else if is_nan_b b then lsbp_qnan (lin b)
  else if (is_zero_b a && is_zero_b b) || (is_inf_b a && is_inf_b b) then lsbp_nan w
  else
    let lm := lin_max w in
    let sc := scale w in
    let a_pos := is_positive_b a || is_zero_b a || is_pos_inf_b a in
    let b_pos := is_positive_b b || is_zero_b b || is_pos_inf_b b in
    let res_pos := Bool.eqb a_pos b_pos in
    if is_zero_b b then
      if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
    else if is_zero_b a then lsbp_zero w
    else if is_inf_b a then if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
    else if is_inf_b b then lsbp_zero w
    else
      let ad := if a_pos then idx a - pit_hi w else pit_lo w - idx a in
      let bd := if b_pos then idx b - pit_hi w else pit_lo w - idx b in
      let ap := (ad * lm + lin a) - lm in
      let bp := (bd * lm + lin b) - lm in
      if Z.eqb bp 0 then if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
      else
        let quot := (ap * sc) / bp + lm in
        let max_pos := pit_lo w * lm in
        if Z.leb max_pos quot then
          if res_pos then lsbp_pos_inf w else lsbp_neg_inf w
        else
          let rd := quot / lm in
          let rl := quot mod lm in
          if Z.eqb rd 0 then lsbp_zero w
          else if res_pos
          then mkLsbp (pit_hi w + rd) rl
          else mkLsbp (pit_lo w - rd) rl.

(** ** Square Root *)

Definition isqrt (n : Z) : Z :=
  if Z.leb n 0 then 0
  else
    let fix newton (x y : Z) (fuel : nat) :=
      match fuel with
      | O => x
      | S f => if Z.ltb y x then newton y ((y + n / y) / 2) f else x
      end
    in newton n ((n + 1) / 2) 64%nat.

Definition lsbp_sqrt {w} (v : Lsbp w) : Lsbp w :=
  if is_nan_b v then lsbp_qnan (lin v)
  else if is_zero_b v then lsbp_zero w
  else if is_pos_inf_b v then lsbp_pos_inf w
  else if is_negative_b v || is_neg_inf_b v then lsbp_nan w
  else
    let lm := lin_max w in
    let sc := scale w in
    let d := idx v - pit_hi w in
    let pos := (d * lm + lin v) - lm in
    let root := isqrt (pos * sc) + lm in
    let rd := root / lm in
    let rl := root mod lm in
    mkLsbp (pit_hi w + rd) rl.

(** ** Absolute Value *)

Definition lsbp_abs {w} (v : Lsbp w) : Lsbp w :=
  if is_negative_b v then lsbp_neg v else v.

(** ** Comparison *)

Definition lsbp_compare {w} (a b : Lsbp w) : comparison :=
  if is_nan_b a || is_nan_b b then Eq
  else Qcompare (decode a) (decode b).

Definition lsbp_eq {w} (a b : Lsbp w) : bool :=
  match lsbp_compare a b with Eq => true | _ => false end.

Definition lsbp_lt {w} (a b : Lsbp w) : bool :=
  match lsbp_compare a b with Lt => true | _ => false end.

Definition lsbp_le {w} (a b : Lsbp w) : bool :=
  match lsbp_compare a b with Lt | Eq => true | _ => false end.
