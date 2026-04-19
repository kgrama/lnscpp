# Build pbf_wasm.cpp → pbf_wasm.{js,wasm} via emscripten.
#
# Two paths:
#   make           — uses `emcc` on PATH (install: pacman -S emscripten, or emsdk)
#   make docker    — runs the official emscripten/emsdk container; no host install
#
# Both produce the same artifacts.

EMCC       ?= emcc
DOCKER_IMG ?= emscripten/emsdk:latest
CXXFLAGS   := -O3 -std=c++11
EXPORTS    := \
  _pbf_init_w _pbf_taper_init_w \
  _pbf_encode_w _pbf_decode_w _pbf_add_w _pbf_sub_w _pbf_mul_w _pbf_div_w \
  _pbf_t_encode_w _pbf_t_decode_w _pbf_t_add_w _pbf_t_sub_w _pbf_t_mul_w _pbf_t_div_w \
  _render_mandelbrot _render_julia _render_burning \
  _malloc _free
RUNTIME    := HEAPU8 HEAPU32 HEAPF64

# turn `_a _b _c` into `["_a","_b","_c"]`
EXPORT_LIST  := $(shell echo '$(EXPORTS)' | sed 's/[^ ][^ ]*/"&"/g; s/ /,/g')
RUNTIME_LIST := $(shell echo '$(RUNTIME)' | sed 's/[^ ][^ ]*/"&"/g; s/ /,/g')

EMCC_CMD = $(EMCC) $(CXXFLAGS) pbf_wasm.cpp -o pbf_wasm.js \
	  -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME=createPbfModule \
	  -s ENVIRONMENT=web -s SINGLE_FILE=0 \
	  -s ALLOW_MEMORY_GROWTH=1 \
	  -s "EXPORTED_FUNCTIONS=[$(EXPORT_LIST)]" \
	  -s "EXPORTED_RUNTIME_METHODS=[$(RUNTIME_LIST)]"

.PHONY: all docker clean serve

all: pbf_wasm.js

# Native build — needs `emcc` on PATH.
pbf_wasm.js: pbf_wasm.cpp ../pbf.cpp ../pbf_taper.cpp
	$(EMCC_CMD)

# Containerised build — needs Docker. Mounts repo root read-write at /src and
# runs emcc with the same flags. First run pulls the image (~600 MB).
docker: pbf_wasm.cpp ../pbf.cpp ../pbf_taper.cpp
	docker run --rm \
	  -v "$$(pwd)/..":/src -w /src/wasm \
	  -u $$(id -u):$$(id -g) \
	  $(DOCKER_IMG) $(EMCC_CMD)

# Quick local server so you can open the demo in a browser.
serve:
	@echo "Open http://localhost:8765/fractal_wasm.html"
	python3 -m http.server 8765

# Containerised static server (nginx:alpine, ~10 MB).
# Foreground: Ctrl-C to stop.
serve-docker:
	@echo "Open http://localhost:8765/fractal_wasm.html"
	docker run --rm -p 8765:80 \
	  -v "$$(pwd)":/usr/share/nginx/html:ro \
	  nginx:alpine

clean:
	rm -f pbf_wasm.js pbf_wasm.wasm
