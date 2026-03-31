SVARP_VERSION := "1.1.1"
SVARP_UPDATE  := "Dec 09, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Build mode
#   USE_CONDA = 0  (default) Build all dependencies from source under dep/.
#                   Run 'make libs' first, then 'make'.
#   USE_CONDA = 1  Use htslib and tools from the conda environment.
#                   Only WFA2 is built locally; 'make' handles it automatically.
# ================================================================
USE_CONDA ?= 0

# Vars
TARGET_EXEC := svarp
BUILD_DIR   := build
SRC_DIRS    := src
DEP_DIR     := dep

HTSLIB_DIR   := $(DEP_DIR)/htslib
WFA_DIR      := $(DEP_DIR)/wfa
WFA_LIB      := $(WFA_DIR)/lib/libwfacpp.a

WTDBG2_DIR   := third_party/wtdbg2
WTDBG2_BIN   := $(WTDBG2_DIR)/wtdbg2

# minimap2 (from source: built under dep/; conda: provided by the environment)
MINIMAP2_DIR := $(DEP_DIR)/minimap2
MINIMAP2_BIN := $(MINIMAP2_DIR)/minimap2

HTSLIB_VERSION := 1.17
HTSLIB_TARBALL := $(DEP_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_LIB     := $(HTSLIB_DIR)/libhts.a

WFA_VERSION     := 2.3.4
WFA_TARBALL     := $(DEP_DIR)/WFA2-lib-$(WFA_VERSION).tar.gz

# ------------------------------------------------
# Compiler / flags
# ------------------------------------------------
CXX ?= g++
BUILD ?= release

CXXFLAGS ?=
CXXFLAGS += -Wall -std=c++17 \
           -DSVARP_VERSION=\"$(SVARP_VERSION)\" \
           -DBUILD_DATE=\"$(BUILD_DATE)\" \
           -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" \
           -DSVARP_DEBUG=$(SVARP_DEBUG)

ifeq ($(BUILD),debug)
    CXXFLAGS += -O0 -g
else
    CXXFLAGS += -O3 -DNDEBUG
endif

LDFLAGS ?=

# ------------------------------------------------
# Include / link paths
# ------------------------------------------------
ifeq ($(USE_CONDA),1)
    CXXFLAGS += -I$(PREFIX)/include -I$(PREFIX)/include/htslib -I$(DEP_DIR)/wfa
    LDFLAGS  += -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib -lhts -lz -lpthread $(WFA_LIB)
else
    CXXFLAGS += -I$(HTSLIB_DIR) -I$(DEP_DIR)/wfa
    LDFLAGS  += $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread
endif

# Sources and objects
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))


.PHONY: all clean clean-libs libs post-build
.PHONY: test-all test

all: $(BUILD_DIR)/$(TARGET_EXEC) post-build

# ================ Tests
test-logfile: build/test_logfile
	./build/test_logfile

build/test_logfile: src/logfile.cpp tests/logfile_test.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I. -Isrc -Idep/htslib -Idep/wfa src/logfile.cpp tests/logfile_test.cpp -o $@

test-common: build/test_common_parse_gaf build/test_cigar build/test_common_utils build/test_run_and_log build/test_variant build/test_generate_sv_node build/test_merge_neighbor_nodes build/test_assembly build/test_remap build/test_alignment build/test_phasing build/test_variant_mapping build/test_alignment_read_gz
	./build/test_common_parse_gaf
	./build/test_cigar
	./build/test_common_utils
	./build/test_run_and_log
	./build/test_variant
	./build/test_generate_sv_node
	./build/test_merge_neighbor_nodes
	./build/test_assembly
	./build/test_remap
	./build/test_alignment
	./build/test_phasing
	./build/test_variant_mapping
	./build/test_alignment_read_gz

test-all: test-logfile test-common
	@echo "=== All tests finished ==="

test: test-all
	@echo "Ran make test (alias to test-all)"

build/test_common_parse_gaf: tests/common_parse_gaf_test.cpp src/common.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/common.cpp src/logfile.cpp tests/common_parse_gaf_test.cpp -o $@

build/test_cigar: tests/cigar_test.cpp src/common.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/common.cpp src/logfile.cpp tests/cigar_test.cpp -o $@

build/test_common_utils: tests/common_utils_test.cpp src/common.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/common.cpp src/logfile.cpp tests/common_utils_test.cpp -o $@

build/test_run_and_log: tests/run_and_log_test.cpp src/common.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/common.cpp src/logfile.cpp tests/run_and_log_test.cpp -o $@


build/test_variant: tests/variant_test.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/variant_test.cpp -o $@


build/test_generate_sv_node: tests/generate_sv_node_test.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/generate_sv_node_test.cpp -o $@


build/test_variant_mapping: tests/variant_mapping_test.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/variant_mapping_test.cpp -o $@


build/test_merge_neighbor_nodes: tests/merge_neighbor_nodes_test.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/merge_neighbor_nodes_test.cpp -o $@


build/test_assembly: tests/assembly_test.cpp src/assembly.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/assembly.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/assembly_test.cpp -o $@ $(LDFLAGS)


build/test_remap: tests/remap_test.cpp src/remap.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/remap.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/remap_test.cpp -o $@ $(LDFLAGS)


build/test_alignment: tests/alignment_test.cpp src/alignment.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/alignment.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/alignment_test.cpp -o $@ $(LDFLAGS)


build/test_phasing: tests/phasing_test.cpp src/phasing.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/phasing.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/phasing_test.cpp -o $@ $(LDFLAGS)


build/test_alignment_read_gz: tests/alignment_read_gz_test.cpp src/alignment.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc -Idep/htslib -Idep/wfa src/alignment.cpp src/variant.cpp src/common.cpp src/reference.cpp src/logfile.cpp tests/alignment_read_gz_test.cpp -o $@ $(LDFLAGS)


# Main binary
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Object compilation rule
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# In conda mode, libs are built automatically before objects.
# Otherwise, run 'make libs' manually first.
ifeq ($(USE_CONDA),1)
$(OBJS): | libs
endif

# ================================================================
# Clean targets
# ================================================================
clean:
	rm -rf $(BUILD_DIR)

clean-libs:
	rm -rf $(DEP_DIR)

# ================================================================
# Build finished message
# ================================================================
post-build:
	@echo ""
	@echo "======================================================"
	@echo " SVarp build finished"
	@echo ""
	@echo " Add to PATH (recommended):"
	@echo ""
	@echo "   export PATH=$(PWD)/build:\$$PATH"
	@echo ""
	@echo " Then you can run:"
	@echo "   svarp"
	@echo "======================================================"
	@echo ""

# ================================================================
# libs target
# ================================================================
ifeq ($(USE_CONDA),1)
# Conda: WFA2 and wtdbg2 built locally; htslib and minimap2 from conda
libs: $(WFA_LIB) $(WTDBG2_BIN)
	@echo "Using conda htslib; WFA2 and wtdbg2 built locally."
	@echo "NOTE: minimap2 and samtools must be available in the conda environment."
else
# From source: all dependencies built locally
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN) $(MINIMAP2_BIN)
	@echo "Built HTSLIB, WFA2, wtdbg2 and minimap2 under dep/"
endif

# ================================================================
# HTSLIB
# ================================================================
$(HTSLIB_LIB): $(HTSLIB_TARBALL)
	mkdir -p $(HTSLIB_DIR)
	tar -xvf $(HTSLIB_TARBALL) -C $(HTSLIB_DIR) --strip-components=1
	cd $(HTSLIB_DIR) && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make

$(HTSLIB_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2 -O $(HTSLIB_TARBALL)

# ================================================================
# WFA2-lib (built under dep/ in both modes)
# ================================================================
$(WFA_LIB): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
ifeq ($(shell uname -s),Darwin)
	@echo "Patching WFA2 endian.h for macOS..."
	find $(WFA_DIR) -name '*.c' -o -name '*.h' | xargs sed -i '' 's|#include <endian.h>|#ifdef __APPLE__\n#include <machine/endian.h>\n#else\n#include <endian.h>\n#endif|g'
endif
	cd $(WFA_DIR) && BUILD_TOOLS=0 BUILD_EXAMPLES=0 $(MAKE) -j1 clean setup lib_wfa

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ================================================================
# minimap2 (skipped if already in PATH)
# ================================================================
$(MINIMAP2_BIN):
	mkdir -p $(DEP_DIR)
	@if command -v minimap2 >/dev/null 2>&1; then \
	  echo "minimap2 found in PATH, skipping local build in dep/minimap2"; \
	else \
	  echo "minimap2 not found in PATH, building under dep/minimap2"; \
	  test -d $(MINIMAP2_DIR) || git clone https://github.com/lh3/minimap2 $(MINIMAP2_DIR); \
	  cd $(MINIMAP2_DIR) && $(MAKE); \
	fi

# ================================================================
# wtdbg2 (bundled in third_party/, requires gcc on macOS)
# ================================================================
WTDBG2_CC := $(CC)
ifeq ($(shell uname -s),Darwin)
  WTDBG2_CC := $(shell command -v gcc-15 2>/dev/null || command -v gcc-14 2>/dev/null || command -v gcc-13 2>/dev/null || echo gcc)
endif

$(WTDBG2_BIN):
	cd $(WTDBG2_DIR) && $(MAKE) CC=$(WTDBG2_CC)
