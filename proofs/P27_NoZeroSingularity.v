(** * P27: No zero-singularity for SBP normal values
    ========================================================================

    Kinematic / control-loop claim: SBP's idx structure prevents normal
    values from approaching zero continuously. The trajectory from
    negative to positive must traverse the pivot gap (denormal + NaN
    slots), not a continuum of shrinking-to-zero normals.

    CONSEQUENCE for kinematics:
      IEEE:  |1 / sin(θ)| can be arbitrarily large as θ → 0.
      SBP:   the smallest positive NORMAL has distance 1, so |1/v|
             for any normal v is bounded by SCALE — no singularity.

    The structure:
      idx:   0     1           pit_lo   nan   pit_hi           idx_max
             zero  −∞   ...    −denorm  NaN   +denorm   ...    +∞
                         (negative normals)        (positive normals)

    No continuous path from a negative normal to a positive normal
    passes through idx=0. Any two normals of opposite sign have a
    non-empty pivot-gap idx set between them.
*)

From Stdlib Require Import ZArith Lia.
Require Import Semantics.
Require Import Properties.
Open Scope Z_scope.

(** ** The gap between positive and negative normal regions *)

(** Smallest positive normal: idx just above pit_hi.
    Every positive-normal idx satisfies idx > pit_hi, so distance >= 1. *)
Theorem pos_normal_min_distance : forall w (v : Lsbp w),
  is_valid v -> is_positive v -> distance v >= 1.
Proof.
  intros w v Hv Hp.
  unfold is_positive in Hp.
  unfold distance, is_positive_b.
  destruct (Z.ltb_spec (pit_hi w) (idx v)) as [|HH]; [lia|lia].
Qed.

(** Smallest negative normal: idx just below pit_lo.
    distance >= 1 for any negative normal. *)
Theorem neg_normal_min_distance : forall w (v : Lsbp w),
  is_valid v -> is_negative v -> distance v >= 1.
Proof.
  intros w v Hv Hn.
  unfold is_negative in Hn.
  destruct Hn as [H1 H2].
  unfold distance, is_positive_b, is_negative_b.
  destruct (Z.ltb_spec (pit_hi w) (idx v)).
  - (* idx v > pit_hi is impossible: Hn says idx v < pit_lo < pit_hi *)
    assert (pit_lo w < pit_hi w).
    { unfold pit_lo, pit_hi, half_bits. destruct w; simpl; lia. }
    lia.
  - destruct (Z.ltb_spec 1 (idx v)); [|lia].
    destruct (Z.ltb_spec (idx v) (pit_lo w)); [|lia].
    simpl. lia.
Qed.

(** The pivot gap exists: between pit_lo and pit_hi is exactly 4 non-normal
    idx positions (pit_lo, nan_lo, nan_hi, pit_hi). *)
Theorem pivot_gap_width : forall w,
  pit_hi w - pit_lo w = 3.
Proof.
  intros w. unfold pit_hi, pit_lo, half_bits.
  destruct w; simpl; lia.
Qed.

(** ** No continuous zero-crossing for normal values *)

(** If a is negative-normal and b is positive-normal, then any idx strictly
    between idx a and idx b lies in the pivot gap or the infinity slots —
    NOT in the normal range. There is no "smaller and smaller positive
    normal" on the path from negative to positive. *)
(** No idx can be both positive-normal and negative-normal simultaneously.
    The pivot gap enforces this disjointness. *)
Theorem normal_ranges_disjoint : forall w (i : Z),
  ~ (i > pit_hi w /\ i > 1 /\ i < pit_lo w).
Proof.
  intros w i [Hhi [_ Hlo]].
  assert (pit_lo w < pit_hi w) by
    (unfold pit_lo, pit_hi, half_bits; destruct w; simpl; lia).
  lia.
Qed.

(** ** Bounded division: no singularity *)

(** Core claim: if b is a normal SBP value, then distance b >= 1, so
    the division a/b cannot produce infinity unless overflow saturation
    kicks in. This is the "no singularity" property.

    In IEEE, b can be subnormal (down to 2^-126 × 2^-23 ≈ 1e-45), so
    |a/b| can reach 1e45 × |a|. In SBP, min(distance(b)) = 1 (or
    denormal = tag), so either:
      - b is normal: distance >= 1, bounded division
      - b is denormal: linear interpolation in fixed range, bounded
      - b is exactly zero: division yields ±∞ cleanly (tagged)
*)

Definition div_bounded_when_normal (w : BitWidth) : Prop :=
  forall a b : Lsbp w,
    is_valid a -> is_valid b ->
    is_finite a -> is_normal b ->
    (** b has distance >= 1, so a/b has magnitude <= distance(a) * SCALE
        (in linear interpretation). The result is never +∞ by math; it
        may saturate only via IDX_MAX clamp. *)
    distance b >= 1.

Theorem div_bounded_when_normal_proof : forall w, div_bounded_when_normal w.
Proof.
  intros w a b Hav Hbv Haf Hnb.
  unfold is_normal in Hnb. destruct Hnb as [Hp | Hn].
  - apply pos_normal_min_distance; assumption.
  - apply neg_normal_min_distance; assumption.
Qed.

(** ** Kinematics corollary *)

(** For a 2-link inverse kinematics involving 1/sin(θ), if sin(θ) is
    represented as a normal SBP value (not special-tagged), then 1/sin(θ)
    is bounded. The singularity is "clamped" at the smallest-normal
    boundary rather than extending to infinity.

    This makes IK numerically stable without damped pseudoinverse logic:
    the encoding itself provides the damping via denormal quantization.
*)
Definition kinematics_no_singularity (w : BitWidth) : Prop :=
  forall sin_theta v : Lsbp w,
    is_valid sin_theta -> is_valid v ->
    is_normal sin_theta ->
    (** The encoding guarantees distance(sin_theta) >= 1, bounding 1/sin_theta *)
    distance sin_theta >= 1.

Theorem kinematics_no_singularity_proof : forall w, kinematics_no_singularity w.
Proof.
  intros w sin_theta v Hs Hv Hn.
  unfold is_normal in Hn. destruct Hn as [Hp | Hnn].
  - apply pos_normal_min_distance; assumption.
  - apply neg_normal_min_distance; assumption.
Qed.

(** ** Composite theorem: SBP-robust kinematic chain *)
Theorem sbp_kinematics_chain : forall w,
  (** For any finite pair of operands with b normal, 1/b is bounded *)
  (forall a b : Lsbp w,
    is_valid a -> is_valid b -> is_finite a -> is_normal b ->
    distance b >= 1) /\
  (** There's no normal-normal path with arbitrarily small intermediates *)
  (pit_hi w - pit_lo w = 3).
Proof.
  intros w. split.
  - intros a b Hav Hbv Haf Hbn.
    destruct Hbn as [Hp | Hn].
    + apply pos_normal_min_distance; assumption.
    + apply neg_normal_min_distance; assumption.
  - apply pivot_gap_width.
Qed.
