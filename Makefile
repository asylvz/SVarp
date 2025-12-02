SVARP_VERSION := "1.0.1"
SVARP_UPDATE  := "Dec 02, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Compiler / build type
# ================================================================
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
    CXXFLAGS += -O3 -DNDEBUG
endif

# ================================================================
# Dependencies (always bare-metal under dep/)
# ================================================================
DEP_DIR       := dep

# HTSLIB
HTSLIB_VERSION := 1.17
HTSLIB_TARBALL := $(DEP_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_DIR     := $(DEP_DIR)/htslib
HTSLIB_LIB     := $(HTSLIB_DIR)/libhts.a

# WFA2-lib
WFA_VERSION    := 2.3.4
WFA_TARBALL    := $(DEP_DIR)/WFA2-lib-$(WFA_VERSION).tar.gz
WFA_DIR        := $(DEP_DIR)/wfa
WFA_LIB        := $(WFA_DIR)/lib/libwfacpp.a
WFA_HEADERS    := $(WFA_DIR)/bindings/cpp/WFAligner.hpp

# wtdbg2
WTDBG2_DIR     := $(DEP_DIR)/wtdbg2
WTDBG2_BIN     := $(WTDBG2_DIR)/wtdbg2

# ================================================================
# Paths / flags
# ================================================================
# Include paths: local htslib + local WFA2
CXXFLAGS += -I$(HTSLIB_DIR) -I$(WFA_DIR)

# Link with local static libs
LDFLAGS  += $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread

# Build directory / sources
TARGET_EXEC := svarp
BUILD_DIR   := build
SRC_DIRS    := src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

.PHONY: all clean clean-libs libs post-build

# ================================================================
# Default target
# ================================================================
all: $(BUILD_DIR)/$(TARGET_EXEC) post-build

# Link step
$(BUILD_DIR)/$(TARGET_EXEC): libs $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# C++ compile rule
# ÖNEMLİ: WFA_HEADERS dependency'si burada; önce WFA2 indir + derle + header'ları aç
$(BUILD_DIR)/%.o: %.cpp $(WFA_HEADERS)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ================================================================
# Libraries (always local under dep/)
# ================================================================
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN)

# ---------------- HTSLIB ----------------
$(HTSLIB_LIB): $(HTSLIB_TARBALL)
	mkdir -p $(HTSLIB_DIR)
	tar -xvf $(HTSLIB_TARBALL) -C $(HTSLIB_DIR) --strip-components=1
	cd $(HTSLIB_DIR) && ./configure --disable-lzma --disable-bz2 --disable-libcurl && $(MAKE)

$(HTSLIB_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2 -O $(HTSLIB_TARBALL)

# ---------------- WFA2-lib ----------------
# Header dosyasını sentinel gibi kullanıyoruz:
$(WFA_HEADERS): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
	cd $(WFA_DIR) && $(MAKE) clean all

$(WFA_LIB): $(WFA_HEADERS)
	@# WFA2 build sonrası lib/libwfacpp.a oluşmuş olmalı
	test -f $(WFA_LIB)

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ---------------- wtdbg2 ----------------
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
	cd $(WTDBG2_DIR) && $(MAKE)

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

