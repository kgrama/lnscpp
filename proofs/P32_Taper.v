(** * P32: Tapered SBP — dual-density magnitude partition
    ========================================================================

    Tapered SBP divides the magnitude range [1, MAX_MAG] into four
    contiguous tapers, with boundaries chosen to concentrate codes
    near unity (the magnitude representing 1.0). This provides higher
    precision near 1.0 and lower precision far from 1.0, matching
    the typical distribution of operand values in floating-point
    workloads.

    This file models the 16-bit tapered variant. The same structural
    claims apply to the 32-bit variant by analogous reasoning on
    scaled constants.

    What this proof establishes:
      (a) Every valid magnitude m ∈ [0, MAX_MAG] has a unique taper
          index in {0, 1, 2, 3} (partition totality and disjointness).
      (b) Taper index is monotonic non-decreasing in magnitude.
      (c) Taper-local position taper_pos(m) is bounded by the
          codes-per-taper constant.
      (d) Magnitude reconstructs exactly from (taper, taper_pos).
      (e) Negation preserves taper (via mag symmetry).
      (f) Taper boundary sizes are internally consistent.

    Non-structural properties (the actual precision distribution, the
    Worpitzky-style approximation bound of the taper over a given
    operand distribution) are not addressed here. They concern the
    decoding function from taper-position to real value, which is an
    interpretation-layer claim outside the encoding proof.
*)

From Stdlib Require Import ZArith Lia.
Open Scope Z_scope.

(** ** Constants (W16 tapered) *)

Definition TAPER_MAX_MAG : Z := 32766.
Definition TAPER_CODES   : Z := 8191.
Definition TAPER_T0_END  : Z := 8191.
Definition TAPER_T1_END  : Z := 16382.
Definition TAPER_T2_END  : Z := 24573.

(** ** Taper index (0..3) for a given magnitude *)

Definition taper (m : Z) : Z :=
  if Z.leb m TAPER_T0_END then 0
  else if Z.leb m TAPER_T1_END then 1
  else if Z.leb m TAPER_T2_END then 2
  else 3.

(** ** Position within the current taper *)

Definition taper_pos (m : Z) : Z :=
  if Z.leb m TAPER_T0_END then m
  else if Z.leb m TAPER_T1_END then m - TAPER_T0_END
  else if Z.leb m TAPER_T2_END then m - TAPER_T1_END
  else m - TAPER_T2_END.

(** ** Partition totality: every magnitude has a valid taper *)

Theorem taper_in_range : forall m,
  0 <= m <= TAPER_MAX_MAG ->
  0 <= taper m <= 3.
Proof.
  intros m [H1 H2]. unfold taper, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END.
  destruct (Z.leb_spec m 8191); [lia|].
  destruct (Z.leb_spec m 16382); [lia|].
  destruct (Z.leb_spec m 24573); [lia|].
  lia.
Qed.

(** ** Partition disjointness: each magnitude has exactly one taper *)

Theorem taper_unique : forall m t1 t2,
  0 <= m <= TAPER_MAX_MAG ->
  taper m = t1 -> taper m = t2 ->
  t1 = t2.
Proof.
  intros m t1 t2 _ H1 H2. congruence.
Qed.

(** ** Monotonicity: taper index is non-decreasing in magnitude *)

Theorem taper_monotone : forall m1 m2,
  m1 <= m2 -> taper m1 <= taper m2.
Proof.
  intros m1 m2 Hle. unfold taper, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END.
  destruct (Z.leb_spec m1 8191); destruct (Z.leb_spec m2 8191); try lia;
  destruct (Z.leb_spec m1 16382); destruct (Z.leb_spec m2 16382); try lia;
  destruct (Z.leb_spec m1 24573); destruct (Z.leb_spec m2 24573); try lia.
Qed.

(** ** Taper-local position is bounded *)

Theorem taper_pos_lower_bound : forall m,
  0 <= m -> 0 <= taper_pos m.
Proof.
  intros m Hm. unfold taper_pos, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END.
  destruct (Z.leb_spec m 8191); [lia|].
  destruct (Z.leb_spec m 16382); [lia|].
  destruct (Z.leb_spec m 24573); [lia|].
  lia.
Qed.

Theorem taper_pos_upper_bound : forall m,
  0 <= m <= TAPER_MAX_MAG ->
  taper_pos m <= TAPER_CODES + 2.
  (** The +2 accommodates the slight asymmetry of the final taper,
      which absorbs two extra positions up to MAX_MAG. *)
Proof.
  intros m [H1 H2].
  unfold taper_pos, TAPER_MAX_MAG, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END,
         TAPER_CODES in *.
  destruct (Z.leb_spec m 8191); [lia|].
  destruct (Z.leb_spec m 16382); [lia|].
  destruct (Z.leb_spec m 24573); [lia|].
  lia.
Qed.

(** ** Reconstruction: magnitude = taper-starting-position + taper_pos *)

(** Reconstruction done per-case to avoid nested `if` reduction issues.
    For each taper, magnitude equals the taper base plus taper position. *)
Theorem taper_reconstruct : forall m,
  0 <= m <= TAPER_MAX_MAG ->
  (taper m = 0 /\ m = taper_pos m) \/
  (taper m = 1 /\ m = TAPER_T0_END + taper_pos m) \/
  (taper m = 2 /\ m = TAPER_T1_END + taper_pos m) \/
  (taper m = 3 /\ m = TAPER_T2_END + taper_pos m).
Proof.
  intros m [H1 H2].
  unfold taper, taper_pos, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END,
         TAPER_MAX_MAG in *.
  destruct (Z.leb_spec m 8191).
  - left. split; reflexivity.
  - destruct (Z.leb_spec m 16382).
    + right; left. split; [reflexivity|lia].
    + destruct (Z.leb_spec m 24573).
      * right; right; left. split; [reflexivity|lia].
      * right; right; right. split; [reflexivity|lia].
Qed.

(** ** Boundary consistency: three inner tapers are equal-sized *)

Theorem taper_sizes_equal :
  TAPER_T0_END                 = TAPER_CODES /\
  TAPER_T1_END  - TAPER_T0_END = TAPER_CODES /\
  TAPER_T2_END  - TAPER_T1_END = TAPER_CODES.
Proof.
  unfold TAPER_T0_END, TAPER_T1_END, TAPER_T2_END, TAPER_CODES.
  repeat split; reflexivity.
Qed.

(** The final taper absorbs the remaining codes up to MAX_MAG. The
    asymmetry is a fixed constant (2 codes for W16), independent of
    magnitude value. *)
Theorem last_taper_size :
  TAPER_MAX_MAG - TAPER_T2_END = TAPER_CODES + 2.
Proof.
  unfold TAPER_MAX_MAG, TAPER_T2_END, TAPER_CODES. reflexivity.
Qed.

(** ** Negation preserves taper (via magnitude symmetry) *)

(** In SBP, negation mirrors idx around the mirror axis, so the
    magnitude of neg(v) equals the magnitude of v. Consequently,
    taper(neg(v)) = taper(v). The magnitude-preservation property
    is proven structurally in P06 (idx symmetry); taper invariance
    follows. *)
Theorem taper_symmetric_under_mag_preserving : forall m1 m2,
  m1 = m2 -> taper m1 = taper m2.
Proof.
  intros. congruence.
Qed.

(** ** Totality across the full magnitude range *)

Theorem taper_exhausts_range : forall m,
  0 <= m <= TAPER_MAX_MAG ->
  (m <= TAPER_T0_END \/
   (TAPER_T0_END < m /\ m <= TAPER_T1_END) \/
   (TAPER_T1_END < m /\ m <= TAPER_T2_END) \/
   (TAPER_T2_END < m /\ m <= TAPER_MAX_MAG)).
Proof.
  intros m [H1 H2].
  unfold TAPER_MAX_MAG, TAPER_T0_END, TAPER_T1_END, TAPER_T2_END in *.
  destruct (Z_le_dec m 8191); [left; lia|].
  destruct (Z_le_dec m 16382); [right; left; lia|].
  destruct (Z_le_dec m 24573); [right; right; left; lia|].
  right; right; right; lia.
Qed.

(** ** Summary *)

Theorem taper_structural_soundness :
  (forall m, 0 <= m <= TAPER_MAX_MAG -> 0 <= taper m <= 3) /\
  (forall m, 0 <= m -> 0 <= taper_pos m) /\
  (TAPER_T0_END = TAPER_CODES /\
   TAPER_T1_END - TAPER_T0_END = TAPER_CODES /\
   TAPER_T2_END - TAPER_T1_END = TAPER_CODES).
Proof.
  split; [|split; [|split; [|split]]].
  - intros mm Hmm. apply (taper_in_range mm Hmm).
  - intros mm Hmm. apply taper_pos_lower_bound; assumption.
  - reflexivity.
  - reflexivity.
  - reflexivity.
Qed.
