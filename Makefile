SVARP_VERSION := "1.0.1"
SVARP_UPDATE  := "Dec 02, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Mode switch (bare-metal vs conda-build)
#   USE_CONDA = 0  -> htslib + WFA2 + wtdbg2 local derlenir
#   USE_CONDA = 1  -> htslib conda'dan (-lhts), WFA2 + wtdbg2 local derlenir
# ================================================================
USE_CONDA ?= 0

# ================================================================
# PATH'e build/ ekle (sadece make süresince)
# ================================================================
export PATH := $(PWD)/build:$(PATH)

# ================================================================
# Dependency directories
# ================================================================
DEP_DIR      := dep

# --- WFA2-lib ---
WFA_DIR      := $(DEP_DIR)/wfa
WFA_VERSION  := 2.3.4
WFA_TARBALL  := $(DEP_DIR)/WFA2-lib-$(WFA_VERSION).tar.gz
WFA_LIB      := $(WFA_DIR)/lib/libwfacpp.a

# --- wtdbg2 ---
WTDBG2_DIR   := $(DEP_DIR)/wtdbg2
WTDBG2_BIN   := $(WTDBG2_DIR)/wtdbg2

# --- htslib (sadece bare-metal modda kullanacağız) ---
HTSLIB_DIR      := $(DEP_DIR)/htslib
HTSLIB_VERSION  := 1.17
HTSLIB_TARBALL  := $(DEP_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_LIB      := $(HTSLIB_DIR)/libhts.a

# ================================================================
# Compiler / linker flags
# ================================================================
CXX ?= g++

BUILD ?= release

CXXFLAGS ?=
CXXFLAGS += -Wall -std=c++17 \
           -DSVARP_VERSION=\"$(SVARP_VERSION)\" \
           -DBUILD_DATE=\"$(BUILD_DATE)\" \
           -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" \
           -DSVAPR_DEBUG=$(SVARP_DEBUG)

# Optimizasyon seviyeleri (build tipine göre)
ifeq ($(BUILD),debug)
    CXXFLAGS += -O0 -g
else
    # Conda ortamında kendi optimizasyon bayraklarını kullanmasına izin ver
    ifneq ($(USE_CONDA),1)
        CXXFLAGS += -O3 -DNDEBUG
    endif
endif

# ------------------------------------------------
# Conda modu: htslib conda'dan, WFA2 + wtdbg2 local
# ------------------------------------------------
ifeq ($(USE_CONDA),1)

    # htslib headerları ve WFA2 için include
    CXXFLAGS += -I$(PREFIX)/include -I$(PREFIX)/include/htslib -I$(WFA_DIR)

    # conda'nın htslib'ine link, WFA2'yi local lib'den al
    LDFLAGS  += -L$(PREFIX)/lib -lhts -lz -lpthread $(WFA_LIB)

else
# ------------------------------------------------
# Bare-metal: htslib + WFA2 + wtdbg2 hepsi local
# ------------------------------------------------
    CXXFLAGS += -I$(HTSLIB_DIR) -I$(WFA_DIR)
    LDFLAGS  += $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread
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

all: $(BUILD_DIR)/$(TARGET_EXEC) post-build

$(BUILD_DIR)/$(TARGET_EXEC): libs $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

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
# Conda'da: htslib'i conda paketi sağlar, biz sadece WFA2 + wtdbg2 derliyoruz
libs: $(WFA_LIB) $(WTDBG2_BIN)
	@echo "Using conda htslib; built local WFA2 and wtdbg2."
else
# Bare metal: htslib + WFA2 + wtdbg2 hepsi local
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
# WFA2-lib (always local)
# ================================================================
$(WFA_LIB): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
	cd $(WFA_DIR) && make clean all

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ================================================================
# wtdbg2 (always local)
# ================================================================
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
	cd $(WTDBG2_DIR) && make

