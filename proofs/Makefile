COQC ?= coqc
COQFLAGS ?= -R . SBP

# Foundation
FOUNDATION := Semantics.v Arithmetic.v Properties.v Tactics.v Layouts.v Interpretations.v

# Structural proofs — properties of the SBP encoding itself
# (independent of arithmetic semantics and interpretation)
PROOFS := P01_NegInvolution.v \
          P06_IdxSymmetry.v \
          P09_SignByPosition.v \
          P10_Denormals.v \
          P11_Constants.v \
          P13_PackUnpack.v \
          P16_Ordering.v \
          P17_Sequence.v \
          P18_BitEfficiency.v \
          P19_NoBarrel.v \
          P21_IdxPartition.v \
          P23_NegSpecials.v \
          P25_InterpretationOrthogonality.v \
          P27_NoZeroSingularity.v \
          P28_PureNoSingularity.v \
          P30_LinGeometric.v \
          P31_PolarAdd.v \
          P32_Taper.v

ALL_V := $(FOUNDATION) $(PROOFS)
ALL_VO := $(ALL_V:.v=.vo)

.PHONY: all clean check

all: $(ALL_VO)
	@echo "All structural proofs verified."

# Foundation
Semantics.vo: Semantics.v
	$(COQC) $(COQFLAGS) $<

Arithmetic.vo: Arithmetic.v Semantics.vo
	$(COQC) $(COQFLAGS) $<

Properties.vo: Properties.v Semantics.vo Arithmetic.vo
	$(COQC) $(COQFLAGS) $<

Tactics.vo: Tactics.v
	$(COQC) $(COQFLAGS) $<

Layouts.vo: Layouts.v Semantics.vo Arithmetic.vo
	$(COQC) $(COQFLAGS) $<

Interpretations.vo: Interpretations.v Semantics.vo
	$(COQC) $(COQFLAGS) $<

# Each proof depends on foundation
$(PROOFS:.v=.vo): %.vo: %.v Semantics.vo Arithmetic.vo Properties.vo Tactics.vo
	@echo "Proving $<..."
	$(COQC) $(COQFLAGS) $<

# P18 additionally depends on P13 (pack injectivity from round-trip)
P18_BitEfficiency.vo: P13_PackUnpack.vo

clean:
	rm -f *.vo *.vok *.vos *.glob .*.aux

check: all
