SVARP_VERSION := "1.1.1"
SVARP_UPDATE  := "Dec 09, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Mode switch (bare-metal vs conda-build)
#   USE_CONDA = 0  -> HTSLIB + WFA2 + WTDGB2 + (opsiyonel) minimap2 local build (dep/)
#                     -> Hiçbir dependency otomatik kurulmaz, hepsi 'make libs' ile.
#   USE_CONDA = 1  -> HTSLIB (ve diğer kütüphaneler) conda'dan,
#                     WFA2 local build (dep/), WTDGB2 ise ortamdan (PATH) beklenir
#                     -> 'make' çağrıldığında sadece WFA2 için libs çalışır.
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

WTDBG2_DIR   := $(DEP_DIR)/wtdbg2
WTDBG2_BIN   := $(WTDBG2_DIR)/wtdbg2

# minimap2 (bare-metal'de dep/ altına kuracağız; conda modunda conda'dan gelir)
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
           -DSVAPR_DEBUG=$(SVARP_DEBUG)

ifeq ($(BUILD),debug)
    CXXFLAGS += -O0 -g
else
    ifneq ($(USE_CONDA),1)
        CXXFLAGS += -O3 -DNDEBUG
    endif
endif

LDFLAGS ?=

# ------------------------------------------------
# Include / link flags (conda vs bare-metal)
# ------------------------------------------------
ifeq ($(USE_CONDA),1)
    # conda build ortamı: htslib / zlib / minimap2 vs PREFIX altında
    CXXFLAGS += -I$(PREFIX)/include -I$(PREFIX)/include/htslib -I$(DEP_DIR)/wfa
    LDFLAGS  += -L$(PREFIX)/lib -lhts -lz -lpthread $(WFA_LIB)
else
    # lokal derleme: htslib ve WFA2'yi dep/ altından derle (veya kullanıcı kendisi ayarlasın)
    CXXFLAGS += -I$(HTSLIB_DIR) -I$(DEP_DIR)/wfa
    LDFLAGS  += $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread
endif

# Kaynak ve obje dosyaları
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


# Asıl binary
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Obje derleme kuralı
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------
# Conda build'de: objelerden önce libs çalışsın
# Bare metal'de: libs otomatik çağrılmasın (kullanıcı make libs desin)
# ------------------------------------------------
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
# conda build'de: sadece WFA2 dep/ altına; htslib, wtdbg2, minimap2 conda env'den
libs: $(WFA_LIB)
	@echo "Using conda HTSLIB (libhts); WFA2 built under dep/."
	@echo "NOTE: wtdbg2 and minimap2 must be installed in the conda environment (e.g. conda install -c bioconda wtdbg2 minimap2)."
else
# bare-metal: tüm dependency'ler sadece make libs ile
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN) $(MINIMAP2_BIN)
	@echo "Built HTSLIB, WFA2, wtdbg2 and minimap2 under dep/"
endif

# ================================================================
# HTSLIB (bare-metal only)
# ================================================================
$(HTSLIB_LIB): $(HTSLIB_TARBALL)
	mkdir -p $(HTSLIB_DIR)
	tar -xvf $(HTSLIB_TARBALL) -C $(HTSLIB_DIR) --strip-components=1
	cd $(HTSLIB_DIR) && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make

$(HTSLIB_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2 -O $(HTSLIB_TARBALL)

# ================================================================
# WFA2-lib (her iki modda da dep/ altına)
# ================================================================
$(WFA_LIB): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
	cd $(WFA_DIR) && BUILD_TOOLS=0 BUILD_EXAMPLES=0 $(MAKE) -j1 clean setup lib_wfa

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ================================================================
# minimap2 (bare-metal: dep/ altına, PATH'te yoksa)
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
# wtdbg2 (dep/ altına, conda'da CC/LDFLAGS override ile)
# ================================================================
ifeq ($(USE_CONDA),0)
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
	cd $(WTDBG2_DIR) && $(MAKE)
endif
