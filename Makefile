CXX      = g++
CXXFLAGS = -std=c++11 -O2 -I.
TEST     = tests/sbp_vs_pbf_test.cpp

all: run

pbf_test:    $(TEST); $(CXX) $(CXXFLAGS)            -o $@ $< -lm
sbp_test:    $(TEST); $(CXX) $(CXXFLAGS) -DUSE_SBP  -o $@ $< -lm
bfloat_test: $(TEST); $(CXX) $(CXXFLAGS) -DUSE_BFLOAT -o $@ $< -lm

run: pbf_test sbp_test bfloat_test
	@echo ""; @echo "══ PBF ══════════════════════════════════════════════"; ./pbf_test
	@echo ""; @echo "══ SBP ══════════════════════════════════════════════"; ./sbp_test
	@echo ""; @echo "══ BFLOAT16 ═════════════════════════════════════════"; ./bfloat_test

clean:
	rm -f pbf_test sbp_test bfloat_test
