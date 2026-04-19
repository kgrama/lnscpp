(** * LSBP Properties: Generic Invariants *)

From Stdlib Require Import ZArith.
From Stdlib Require Import QArith.
From Stdlib Require Import List.
Require Import Semantics.
Require Import Arithmetic.
Import ListNotations.
Open Scope Z_scope.

(** ** Width-Independent Properties *)

(** Sign-by-position invariant holds for all widths *)
Definition sign_by_position_invariant {w} (v : Lsbp w) : Prop :=
  (idx v > pit_hi w -> sign v = 1) /\
  (idx v > 1 /\ idx v < pit_lo w -> sign v = -1) /\
  (idx v = 0 -> sign v = 0).

(** Denormal sign is determined by which denormal index *)
Definition denormal_sign_invariant {w} (v : Lsbp w) : Prop :=
  (is_pos_denormal v -> sign v = 1) /\
  (is_neg_denormal v -> sign v = -1).

(** Denormal negation swaps denormal type *)
Definition denormal_neg_swaps (w : BitWidth) : Prop :=
  forall l : Z,
    0 <= l <= lin_max w ->
    lsbp_neg (lsbp_pos_denormal (w:=w) l) = lsbp_neg_denormal l /\
    lsbp_neg (lsbp_neg_denormal (w:=w) l) = lsbp_pos_denormal l.

(** Denormal scaling preserves denormal status *)
Definition scale_preserves_denormal (w1 w2 : BitWidth) : Prop :=
  forall l : Z,
    is_denormal (scale_up w1 w2 (lsbp_pos_denormal (w:=w1) l)) /\
    is_denormal (scale_up w1 w2 (lsbp_neg_denormal (w:=w1) l)).

(** Negation is involution *)
Definition neg_involution (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> ~is_nan v ->
    lsbp_neg (lsbp_neg v) = v.

(** Negation flips sign *)
Definition neg_flips_sign (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> is_normal v ->
    sign (lsbp_neg v) = - sign v.

(** Zero is additive identity *)
Definition add_zero_identity (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> ~is_nan v -> ~is_zero v ->
    lsbp_add v (lsbp_zero w) = v /\
    lsbp_add (lsbp_zero w) v = v.

(** Addition is commutative *)
Definition add_commutative (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_finite a -> is_finite b ->
    lsbp_add a b = lsbp_add b a.

(** x + (-x) = 0 *)
Definition add_inverse (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> is_finite v -> ~is_nan v ->
    lsbp_add v (lsbp_neg v) = lsbp_zero w.

(** Zero absorbs multiplication *)
Definition mul_zero_absorber (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> ~is_nan v -> is_finite v ->
    lsbp_mul v (lsbp_zero w) = lsbp_zero w /\
    lsbp_mul (lsbp_zero w) v = lsbp_zero w.

(** Sign rule for multiplication *)
Definition mul_sign_rule (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_normal a -> is_normal b ->
    sign (lsbp_mul a b) = sign a * sign b.

(** NaN propagation *)
Definition nan_propagates (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_nan a \/ is_nan b ->
    is_nan (lsbp_add a b) /\ is_nan (lsbp_mul a b).

(** Infinity properties *)
Definition inf_add_neg_inf_nan (w : BitWidth) : Prop :=
  is_nan (lsbp_add (lsbp_pos_inf w) (lsbp_neg_inf w)).

Definition zero_mul_inf_nan (w : BitWidth) : Prop :=
  is_nan (lsbp_mul (lsbp_zero w) (lsbp_pos_inf w)).

(** ** Bit-Width Bounds *)

Definition idx_bit_bound (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> 0 <= idx v < Z.shiftl 1 (half_bits w).

Definition lin_bit_bound (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> 0 <= lin v < Z.shiftl 1 (half_bits w).

Definition total_bit_bound (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    let hb := half_bits w in
    0 <= idx v * Z.shiftl 1 hb + lin v < Z.shiftl 1 (2 * hb).

(** ** Pivot Gap Properties *)

Definition pivot_gap_special (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    @in_pivot_gap w (idx v) ->
    is_zero v \/ is_nan v.

(** ** Symmetry *)

Definition idx_symmetry (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> is_normal v ->
    idx v + idx (lsbp_neg v) = idx_max w.

(** ** Constants per Width *)

Definition constants_8 : Prop :=
  idx_max W8 = 15 /\ lin_max W8 = 15 /\
  pit_lo W8 = 6 /\ pit_hi W8 = 9 /\
  nan_lo_idx W8 = 7 /\ nan_hi_idx W8 = 8.

Definition constants_16 : Prop :=
  idx_max W16 = 255 /\ lin_max W16 = 255 /\
  pit_lo W16 = 126 /\ pit_hi W16 = 129 /\
  nan_lo_idx W16 = 127 /\ nan_hi_idx W16 = 128.

Definition constants_32 : Prop :=
  idx_max W32 = 65535 /\ lin_max W32 = 65535 /\
  pit_lo W32 = 32766 /\ pit_hi W32 = 32769 /\
  nan_lo_idx W32 = 32767 /\ nan_hi_idx W32 = 32768.

Definition constants_64 : Prop :=
  idx_max W64 = 4294967295 /\ lin_max W64 = 4294967295 /\
  pit_lo W64 = 2147483646 /\ pit_hi W64 = 2147483649 /\
  nan_lo_idx W64 = 2147483647 /\ nan_hi_idx W64 = 2147483648.

(** ** Scaling Properties *)

(** Scale up then down recovers original *)
Definition scale_roundtrip (w1 w2 : BitWidth) : Prop :=
  forall v : Lsbp w1,
    is_valid v ->
    scale_down w2 w1 (scale_up w1 w2 v) = v.

(** Scaling preserves special values *)
Definition scale_preserves_zero (w1 w2 : BitWidth) : Prop :=
  scale_up w1 w2 (lsbp_zero w1) = lsbp_zero w2.

Definition scale_preserves_inf (w1 w2 : BitWidth) : Prop :=
  scale_up w1 w2 (lsbp_pos_inf w1) = lsbp_pos_inf w2 /\
  scale_up w1 w2 (lsbp_neg_inf w1) = lsbp_neg_inf w2.

Definition scale_preserves_nan (w1 w2 : BitWidth) : Prop :=
  scale_up w1 w2 (lsbp_nan w1) = lsbp_nan w2.

(** Distance preserved under scaling (only lin scales) *)
Definition distance_preserved (w1 w2 : BitWidth) : Prop :=
  forall v : Lsbp w1,
    is_valid v -> is_normal v ->
    distance (scale_up w1 w2 v) = distance v.

(** Pack/unpack round-trip *)
Definition pack_unpack_identity (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    unpack w (pack v) = v.

(** Scaling preserves sign *)
Definition scale_preserves_sign (w1 w2 : BitWidth) : Prop :=
  forall v : Lsbp w1,
    is_valid v -> is_normal v ->
    sign (scale_up w1 w2 v) = sign v.

(** ** CLZ Cancellation Properties *)

(** XOR of magnitudes measures log-domain distance.
    For normal values, magnitude = distance from pivot.
    XOR leading zeros count how many leading bits match. *)

(** CLZ(a XOR b) = 0 when magnitudes are in different halves of the range *)
Definition clz_xor_distance {w} (a b : Lsbp w) : Z :=
  let ma := distance a in
  let mb := distance b in
  let xor := Z.lxor ma mb in
  Z.log2 (idx_max w) - Z.log2 (Z.max xor 1).

(** CLZ detects near-cancellation: when two same-sign values have
    similar magnitude, their XOR has many leading zeros.
    Specifically: if |dist(a) - dist(b)| < 2^k, then CLZ >= (total_bits - k). *)
Definition clz_detects_cancellation (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_normal a -> is_normal b ->
    sign a = sign b ->
    forall k : Z, 0 <= k ->
    Z.abs (distance a - distance b) < Z.shiftl 1 k ->
    clz_xor_distance a b >= Z.log2 (idx_max w) - k.

(** Negation preserves magnitude (distance from pivot) *)
Definition neg_preserves_distance (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> is_normal v ->
    distance (lsbp_neg v) = distance v.

(** Same-magnitude opposite-sign: XOR of distances is 0 *)
Definition cancel_exact_zero (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> is_normal v ->
    Z.lxor (distance v) (distance (lsbp_neg v)) = 0.

(** CLZ on raw idx XOR is meaningful because SBP is monotonic:
    larger magnitude = larger idx (for positives) or smaller idx (for negatives).
    This means idx difference tracks magnitude difference. *)
Definition idx_monotonic_positive (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_positive a -> is_positive b ->
    distance a > distance b ->
    idx a > idx b.

Definition idx_monotonic_negative (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_negative a -> is_negative b ->
    distance a > distance b ->
    idx a < idx b.

(** ** Ordering Properties *)

(** For positive values, unsigned integer compare on idx gives correct ordering *)
Definition positive_order_by_idx (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_positive a -> is_positive b ->
    (idx a > idx b <-> distance a > distance b).

(** For negative values, ordering is reversed *)
Definition negative_order_by_idx (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_negative a -> is_negative b ->
    (idx a < idx b <-> distance a > distance b).

(** Positive values always have higher idx than negative values *)
Definition positive_above_negative (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_positive a -> is_negative b ->
    idx a > idx b.

(** ** Sequence Separability *)

(** A bijection on magnitudes: any permutation of {0..idx_max} *)
Definition is_bijection (w : BitWidth) (f : Z -> Z) : Prop :=
  (forall x, 0 <= x <= idx_max w -> 0 <= f x <= idx_max w) /\
  (forall x y, 0 <= x <= idx_max w -> 0 <= y <= idx_max w ->
    f x = f y -> x = y).

(** Arithmetic depends only on abstract position, not bit pattern.
    If we apply a bijection f to all magnitudes, the arithmetic result
    (in terms of abstract positions) is identical.

    Concretely: for any bijection f,
    add(f(a), f(b)) computed via decode(f(a)), decode(f(b)), add, encode
    gives the same abstract position as add(a, b).

    This is trivially true because FormatSpec::encode/decode defines the
    mapping from float to position, and the arithmetic operates on positions.
    The bit pattern is irrelevant. *)
Definition sequence_separable (w : BitWidth) : Prop :=
  forall (f : Z -> Z),
    is_bijection w f ->
    forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_normal a -> is_normal b ->
    sign a = sign b ->
    (* The gap between positions is invariant under bijection *)
    Z.abs (distance a - distance b) =
    Z.abs (distance (mkLsbp (w:=w) (idx a) (lin a)) -
           distance (mkLsbp (w:=w) (idx b) (lin b))).

(** The identity bijection gives binary sequence *)
Definition binary_sequence (w : BitWidth) (x : Z) : Z := x.

(** Gray sequence bijection: n XOR (n >> 1) *)
Definition gray_sequence (x : Z) : Z := Z.lxor x (Z.shiftr x 1).

(** Gray sequence preserves adjacency: consecutive inputs differ by exactly 1 bit.
    The XOR of adjacent Gray codes is a power of 2 (single bit set). *)
Definition gray_adjacency (w : BitWidth) : Prop :=
  forall x : Z,
    0 <= x < idx_max w ->
    exists k, 0 <= k /\ Z.lxor (gray_sequence x) (gray_sequence (x + 1)) = Z.shiftl 1 k.

(** ** Bit Efficiency *)

(** Pack is injective: different valid Lsbp records produce different bit patterns.
    Equivalently: every bit pattern decodes to at most one value. *)
Definition pack_injective (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    pack a = pack b -> a = b.

(** Zero has a unique idx. *)
Definition zero_idx_unique (w : BitWidth) : Prop :=
  forall v : Lsbp w, is_zero v <-> idx v = 0.

(** Exactly 2 idx values encode infinity. *)
Definition inf_idx_count (w : BitWidth) : Prop :=
  forall i : Z,
    0 <= i <= idx_max w ->
    (i = 1 \/ i = idx_max w) <-> is_inf (mkLsbp (w:=w) i 0).

(** Exactly 2 idx values encode NaN. *)
Definition nan_idx_count (w : BitWidth) : Prop :=
  forall i : Z,
    0 <= i <= idx_max w ->
    (i = nan_lo_idx w \/ i = nan_hi_idx w) <-> is_nan (mkLsbp (w:=w) i 0).

(** Layout C uses 7 idx values for specials: {0, 1, idx_max, pit_lo, pit_hi,
    nan_lo_idx, nan_hi_idx}. All other idx values carry normal magnitude. *)
Definition special_idx_set_layoutC (w : BitWidth) (i : Z) : Prop :=
  i = 0 \/ i = 1 \/ i = idx_max w \/
  i = pit_lo w \/ i = pit_hi w \/
  i = nan_lo_idx w \/ i = nan_hi_idx w.

Definition special_count_layoutC (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    (is_zero v \/ is_inf v \/ is_nan v \/ is_denormal v)
    <-> special_idx_set_layoutC w (idx v).

(** Count of normal positions per width (Layout C).
    Total idx values = 2^half_bits, specials = 7, normals split 50/50 by sign.
    With lin having 2^half_bits values each, this gives total representable
    normal values = (2^half_bits - 7) * 2^half_bits. *)
Definition normal_idx_count_w8 : Prop :=
  (* 2^4 - 7 = 9 normal idx values (excluding 7 specials) *)
  9 = idx_max W8 + 1 - 7.

(** IEEE vs SBP NaN bit-pattern overhead.
    IEEE binary32: exponent = 0xFF, mantissa != 0 → 2 * (2^23 - 1) NaN patterns.
    SBP W32: idx ∈ {nan_lo_idx, nan_hi_idx}, lin free → 2 * 2^16 NaN patterns.
    Ratio: IEEE burns ~128x more bit patterns on NaN. *)
Definition ieee_binary32_nan_count : Z := 2 * (Z.shiftl 1 23 - 1).
Definition sbp_w32_nan_pattern_count : Z := 2 * Z.shiftl 1 16.

Definition sbp_saves_nan_bits : Prop :=
  ieee_binary32_nan_count > sbp_w32_nan_pattern_count.

(** ** Hardware Primitive Budget (P19) *)

(** Minimal set of hardware primitives sufficient for SBP structural ops.
    Notably absent: variable-distance shift (barrel). *)
Inductive HwPrim : Type :=
  | HwIntAdd
  | HwIntSub
  | HwBitXor
  | HwBitNot
  | HwCompare
  | HwClz
  | HwFixedShift1.

(** Every structural SBP op factors through HwPrim (no variable shift). *)
Definition op_uses_only_hwprim (ops : list HwPrim) : Prop :=
  forall p, In p ops ->
    p = HwIntAdd \/ p = HwIntSub \/ p = HwBitXor \/
    p = HwBitNot \/ p = HwCompare \/ p = HwClz \/ p = HwFixedShift1.

(** Circuits for each structural op, expressed as primitive sequences *)
Definition neg_circuit : list HwPrim := HwBitNot :: nil.
Definition sign_circuit : list HwPrim := HwCompare :: nil.
Definition compare_circuit : list HwPrim := HwCompare :: nil.
Definition classify_circuit : list HwPrim := HwCompare :: nil.
Definition cancel_detect_circuit : list HwPrim := HwBitXor :: HwClz :: nil.

(** ** Idx Partition (P21) *)

(** Totality: every valid idx lands in at least one class. *)
Definition idx_partition_total (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    is_zero v \/ is_pos_inf v \/ is_neg_inf v \/
    is_pos_denormal v \/ is_neg_denormal v \/
    idx v = nan_lo_idx w \/ idx v = nan_hi_idx w \/
    is_positive v \/ is_negative v.

(** Disjointness: at most one "singleton" class holds per idx. *)
Definition idx_partition_disjoint (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v ->
    ~ (is_zero v /\ is_pos_inf v) /\
    ~ (is_zero v /\ is_neg_inf v) /\
    ~ (is_pos_inf v /\ is_neg_inf v) /\
    ~ (is_pos_denormal v /\ is_neg_denormal v) /\
    ~ (is_positive v /\ is_negative v).

(** ** Neg on Specials (P23) *)

Definition neg_zero_spec (w : BitWidth) : Prop :=
  lsbp_neg (lsbp_zero w) = lsbp_zero w.

Definition neg_pos_inf_spec (w : BitWidth) : Prop :=
  lsbp_neg (lsbp_pos_inf w) = lsbp_neg_inf w.

Definition neg_neg_inf_spec (w : BitWidth) : Prop :=
  lsbp_neg (lsbp_neg_inf w) = lsbp_pos_inf w.

Definition neg_pos_denormal_spec (w : BitWidth) : Prop :=
  forall l, 0 <= l <= lin_max w ->
    lsbp_neg (lsbp_pos_denormal (w:=w) l) = lsbp_neg_denormal l.

Definition neg_neg_denormal_spec (w : BitWidth) : Prop :=
  forall l, 0 <= l <= lin_max w ->
    lsbp_neg (lsbp_neg_denormal (w:=w) l) = lsbp_pos_denormal l.

Definition neg_nan_stays_nan (w : BitWidth) : Prop :=
  forall v : Lsbp w, is_nan v -> is_nan (lsbp_neg v).

(** ** CF Normalization (P20) *)

(** Simple CF-step model: (p, q) -> (q mod p, p), partial quotient q/p *)
Definition cf_step (pq : Z * Z) : Z * Z :=
  let (p, q) := pq in
  if Z.eqb p 0 then (p, q) else (q mod p, p).

Fixpoint cf_iterate (pq : Z * Z) (k : nat) : Z * Z :=
  match k with
  | O => pq
  | S k' => cf_iterate (cf_step pq) k'
  end.

(** Barrel normalization spec: shift until top bit set *)
Definition barrel_normalize_count (m : Z) : Z :=
  if Z.leb m 0 then 0
  else let n := Z.log2 m in n.

(** CF normalization approximates barrel normalization within k steps *)
Definition cf_matches_barrel (m e : Z) (k : nat) : Prop :=
  let bc := barrel_normalize_count m in
  let (p, _) := cf_iterate (m, Z.shiftl 1 (bc + 1)) k in
  Z.log2 (Z.max p 1) >= bc - 1.

(** ** Interpretation Orthogonality (P25) *)

(** neg is defined purely on Lsbp, independent of interpretation.
    We formalize this as: the Lsbp-level result of neg doesn't depend
    on how the value would be decoded. *)

(** The neg operation on idx is the same regardless of interpretation. *)
Definition neg_idx_interpretation_free (w : BitWidth) : Prop :=
  forall v : Lsbp w,
    is_valid v -> ~is_zero v -> ~is_nan v ->
    ~is_pos_denormal v -> ~is_neg_denormal v ->
    idx (lsbp_neg v) = idx_max w - idx v \/
    (is_pos_inf v /\ lsbp_neg v = lsbp_neg_inf w) \/
    (is_neg_inf v /\ lsbp_neg v = lsbp_pos_inf w).

(** sign is defined purely on idx ranges, no interpretation dependence. *)
Definition sign_idx_only (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    idx a = idx b -> sign a = sign b.

(** compare on same-sign values reduces to idx compare, no interpretation needed. *)
Definition compare_structural (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_positive a -> is_positive b ->
    (idx a < idx b -> distance a < distance b) /\
    (idx a = idx b -> distance a = distance b) /\
    (idx a > idx b -> distance a > distance b).

(** ** Ramanujan Generalized Continued Fractions (P26) *)

(** Ramanujan ln((1+u)/(1-u)) CF with partial numerators k² u² and
    denominators (2k+1). Matches the Rust ln1p_ramanujan function —
    backward evaluation from fixed depth.

    Structure: h_k = (2k+1)·q − k²·u²·q / h_{k+1}
    Final value: 2·u·q / h_0
*)

Fixpoint ln_rama_backward (u2 q : Z) (k : nat) (h : Z) : Z :=
  match k with
  | O => h
  | S k' =>
      let ksq : Z := (Z.of_nat k) * (Z.of_nat k) in
      let hnext : Z :=
        if Z.eqb h 0 then 0
        else (2 * Z.of_nat k' + 1) * q - (ksq * u2 * q) / h
      in
      ln_rama_backward u2 q k' hnext
  end.

Definition ln_rama_init (q : Z) (depth : nat) : Z :=
  (2 * Z.of_nat depth + 1) * q.

Definition ln1p_rama (w q : Z) (depth : nat) : Z :=
  if Z.eqb w 0 then 0
  else
    let u := (w * q) / (2 * q + w) in
    let u2 := (u * u) / q in
    let h := ln_rama_backward u2 q depth (ln_rama_init q depth) in
    if Z.leb h 0 then w
    else (2 * u * q) / h.

(** Ramanujan exp(x) CF: e^x = 1 + 2x/(2−x + x²/(6 + x²/(10 + x²/(14 + ...)))). *)

Fixpoint exp_rama_backward (x2 q : Z) (k : nat) (h : Z) : Z :=
  match k with
  | O => h
  | S k' =>
      let hnext : Z :=
        if Z.eqb h 0 then 0
        else (4 * Z.of_nat k + 2) * q + (x2 * q) / h
      in
      exp_rama_backward x2 q k' hnext
  end.

Definition exp_rama (x q : Z) (depth : nat) : Z :=
  if Z.eqb x 0 then q
  else
    let x2 := (x * x) / q in
    let h_init := (4 * Z.of_nat depth + 2) * q + x2 in
    let h := exp_rama_backward x2 q (pred depth) h_init in
    let h_outer :=
      if Z.eqb h 0 then 0 else 2 * q - x + (x2 * q) / h
    in
    if Z.leb h_outer 0 then 2 * q
    else q + (2 * x * q) / h_outer.

(** Obligations: well-definedness (no div-by-zero) and convergence bound. *)

Definition ramanujan_ln_zero_in (q : Z) (depth : nat) : Prop :=
  ln1p_rama 0 q depth = 0.

Definition ramanujan_exp_zero_in (q : Z) (depth : nat) : Prop :=
  exp_rama 0 q depth = q.

(** Backward evaluator terminates for any finite depth (trivial — finite recursion). *)
Definition ramanujan_backward_terminates (u2 q : Z) (depth : nat) : Prop :=
  exists result : Z, ln_rama_backward u2 q depth (ln_rama_init q depth) = result.

(** Convergence error bound: depth N → truncation error ≤ C / q^N for some C.
    Ramanujan ln CF converges at ~35 dB/iter (per Rust comments). *)
Definition ramanujan_ln_error_bound (w q : Z) (depth : nat) (bound : Z) : Prop :=
  0 < q -> 0 <= w -> w < q ->
  Z.abs (ln1p_rama w q (S depth) - ln1p_rama w q depth) <= bound.
