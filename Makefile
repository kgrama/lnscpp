CXX      = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I.
LDLIBS   = -lm
BUILD    = build

# ── PBF (Parameterized Bounded Format) targets ────────────────────────────────
PBF_BINS = $(BUILD)/pbftest $(BUILD)/pbf_xlns_test $(BUILD)/pbf_taper_test $(BUILD)/pbf_bench

# ── xlns32 (32-bit LNS) targets ───────────────────────────────────────────────
XLNS32_BINS = $(BUILD)/xlns32test $(BUILD)/xlns32funtest \
              $(BUILD)/xlns32_new_functions_test \
              $(BUILD)/xlns32_cmp_util_const_test \
              $(BUILD)/xlns32_explog_test \
              $(BUILD)/xlns32_batch_layernorm_test

all: $(PBF_BINS) $(XLNS32_BINS)

$(BUILD):
	@mkdir -p $(BUILD)

# ── PBF builds ────────────────────────────────────────────────────────────────
$(BUILD)/pbftest:       pbftest.cpp pbf.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

$(BUILD)/pbf_xlns_test: pbf_xlns_test.cpp pbf_xlns.cpp pbf_batch.cpp pbf.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

$(BUILD)/pbf_taper_test: pbf_taper_test.cpp pbf_taper.cpp pbf.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

$(BUILD)/pbf_bench: pbf_bench.cpp pbf_taper.cpp pbf.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

# ── xlns32 top-level tests ────────────────────────────────────────────────────
$(BUILD)/xlns32test:                   xlns32test.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
$(BUILD)/xlns32funtest:                xlns32funtest.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
$(BUILD)/xlns32_new_functions_test:    xlns32_new_functions_test.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

# ── xlns32 tests/ subdirectory ────────────────────────────────────────────────
$(BUILD)/xlns32_cmp_util_const_test:   tests/xlns32_cmp_util_const_test.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
$(BUILD)/xlns32_explog_test:           tests/xlns32_explog_test.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
$(BUILD)/xlns32_batch_layernorm_test:  tests/xlns32_batch_layernorm_test.cpp xlns32.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

# ── Convenience targets ───────────────────────────────────────────────────────
.PHONY: pbf xlns32 run clean

pbf:    $(PBF_BINS)
xlns32: $(XLNS32_BINS)

run: $(BUILD)/pbftest $(BUILD)/pbf_xlns_test $(BUILD)/pbf_taper_test $(BUILD)/pbf_bench
	@echo ""; echo "══ PBF SNR demo ═════════════════════════════════════"
	$(BUILD)/pbftest
	@echo ""; echo "══ Tapered PBF SNR (denser near unity) ══════════════"
	$(BUILD)/pbf_taper_test
	@echo ""; echo "══ Speed bench: PBF vs xlns32 ═══════════════════════"
	$(BUILD)/pbf_bench
	@echo ""; echo "══ xlns16 API shim over PBF16 ═══════════════════════"
	$(BUILD)/pbf_xlns_test

clean:
	rm -rf $(BUILD)
