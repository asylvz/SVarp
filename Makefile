SVARP_VERSION := "1.0.1"
SVARP_UPDATE  := "Nov 30, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Mode switch (bare-metal vs conda-build)
#   USE_CONDA = 0  -> HTSLIB + WFA2 + WTDGB2 local build (GitHub)
#   USE_CONDA = 1  -> HTSLIB + WFA2 conda'dan, WTDGB2 conda'dan (Bioconda)
# ================================================================
USE_CONDA ?= 0

# ================================================================
# Add build/ to PATH during make
# (only affects child processes, not user shell)
# ================================================================
export PATH := $(PWD)/build:$(PATH)

# ================================================================
# Dependency directories
# ================================================================
DEP_DIR      := dep
HTSLIB_DIR   := $(DEP_DIR)/htslib
WFA_DIR      := $(DEP_DIR)/wfa
WFA_LIB      := $(WFA_DIR)/lib/libwfacpp.a

WTDBG2_DIR   := $(DEP_DIR)/wtdbg2
WTDBG2_BIN   := $(WTDBG2_DIR)/wtdbg2

HTSLIB_VERSION  := 1.17
HTSLIB_TARBALL  := $(DEP_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_LIB      := $(HTSLIB_DIR)/libhts.a

WFA_VERSION     := 2.3.4
WFA_TARBALL     := $(DEP_DIR)/WFA2-lib-$(WFA_VERSION).tar.gz

# ================================================================
# Compiler / linker flags
# ================================================================
CXX = g++

# Build type: release (default) veya debug
BUILD ?= release

CXXFLAGS = -Wall -std=c++17 \
           -DSVARP_VERSION=\"$(SVARP_VERSION)\" \
           -DBUILD_DATE=\"$(BUILD_DATE)\" \
           -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" \
           -DSVAPR_DEBUG=$(SVARP_DEBUG)

# Optimizasyon seviyeleri (build tipine göre)
ifeq ($(BUILD),debug)
    CXXFLAGS += -O0 -g
else
    CXXFLAGS += -O3 -DNDEBUG
endif

# ------------------------------------------------
# Conda modu (htslib + wfa2-lib conda'dan)
# ------------------------------------------------
ifeq ($(USE_CONDA),1)

    # Header'lar conda env'den
    CXXFLAGS += -I$(PREFIX)/include -I$(PREFIX)/include/wfa2lib

    # WFA2-lib kütüphanesini dosya bazında bul (isim ne olursa olsun):
    WFA2_LIB := $(firstword \
        $(wildcard $(PREFIX)/lib/libwfa*.so) \
        $(wildcard $(PREFIX)/lib/libwfa*.a)  \
    )

    ifeq ($(WFA2_LIB),)
        $(error "Could not find WFA2-lib library under $(PREFIX)/lib (tried libwfa*.so / libwfa*.a)")
    endif

    # Link conda kütüphaneleri (dosya yolu ile)
    LDFLAGS  = -L$(PREFIX)/lib $(WFA2_LIB) -lhts -lz -lpthread

else
# ------------------------------------------------
# Bare-metal modu (htslib + WFA2 + wtdbg2 local)
# ------------------------------------------------
    CXXFLAGS += -I$(HTSLIB_DIR) -I$(WFA_DIR)
    LDFLAGS  = $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread
endif


# ================================================================
# Build targets
# ================================================================
TARGET_EXEC := svarp
BUILD_DIR   := ./build
SRC_DIRS    := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

.PHONY: all clean clean-libs libs post-build

# Main target
all: $(BUILD_DIR)/$(TARGET_EXEC) post-build

$(BUILD_DIR)/$(TARGET_EXEC): libs $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# C++ compile rule
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

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
libs:
	@echo "Using conda-provided libraries (htslib, wfa2-lib, wtdbg); nothing to build."
else
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN)
endif

# ================================================================
# HTSLIB (only bare-metal)
# ================================================================
$(HTSLIB_LIB): $(HTSLIB_TARBALL)
	mkdir -p $(HTSLIB_DIR)
	tar -xvf $(HTSLIB_TARBALL) -C $(HTSLIB_DIR) --strip-components=1
	cd $(HTSLIB_DIR) && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make

$(HTSLIB_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2 -O $(HTSLIB_TARBALL)

# ================================================================
# WFA2-lib (only bare-metal)
# ================================================================
$(WFA_LIB): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
	cd $(WFA_DIR) && make clean all

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ================================================================
# wtdbg2 (bare-metal: dep/ altına git clone)
# ================================================================
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
	cd $(WTDBG2_DIR) && make

