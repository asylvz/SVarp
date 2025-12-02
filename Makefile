SVARP_VERSION := "1.0.1"
SVARP_UPDATE  := "Dec 02, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ================================================================
# Mode switch (bare-metal vs conda-build)
#   USE_CONDA = 0  -> HTSLIB + WFA2 + WTDGB2 local build (dep/)
#   USE_CONDA = 1  -> HTSLIB + WFA2 conda env'den, WTDGB2 dep/ altına conda CC ile
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
    # conda build ortamı: htslib / wfa2-lib / zlib vs zaten PREFIX altında
    CXXFLAGS += -I$(PREFIX)/include -I$(PREFIX)/include/htslib -I$(DEP_DIR)/wfa
    LDFLAGS  += -L$(PREFIX)/lib -lhts -lz -lpthread $(WFA_LIB)
else
    # lokal derleme: htslib ve WFA2'yi dep/ altından derle
    CXXFLAGS += -I$(HTSLIB_DIR) -I$(DEP_DIR)/wfa
    LDFLAGS  += $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread
endif

# Kaynak ve obje dosyaları
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))


.PHONY: all clean clean-libs libs post-build

all: $(BUILD_DIR)/$(TARGET_EXEC) post-build

# Asıl binary
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Obje derleme kuralı
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Önce libs, sonra objeler (race condition fix)
$(OBJS): | libs

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
# conda build'de: WFA2'yi yine dep/ altına indirip statik lib derliyoruz,
# ama htslib'ı conda'dan kullanıyoruz. wtdbg2 her iki modda da dep/ altına.
libs: $(WFA_LIB) $(WTDBG2_BIN)
	@echo "Using conda HTSLIB (libhts) from PREFIX, WFA2 and wtdbg2 from dep/"
else
# bare-metal: htslib + WFA2 + wtdbg2'yi dep/ altına full derle
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN)
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
# wtdbg2 (dep/ altına, conda'da CC/LDFLAGS override ile)
# ================================================================
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
ifeq ($(USE_CONDA),1)
	cd $(WTDBG2_DIR) && $(MAKE) \
	    CC="$(BUILD_PREFIX)/bin/x86_64-conda-linux-gnu-cc" \
	    CFLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem $(PREFIX)/include -fdebug-prefix-map=$(SRC_DIR)=/usr/local/src/conda/svarp-$(SVARP_VERSION) -fdebug-prefix-map=$(PREFIX)=/usr/local/src/conda-prefix" \
	    CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem $(PREFIX)/include" \
	    LDFLAGS="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,$(PREFIX)/lib -Wl,-rpath-link,$(PREFIX)/lib -L$(PREFIX)/lib -L$(PREFIX)/lib -lhts -lz -lpthread $(WFA_LIB)"
else
	cd $(WTDBG2_DIR) && $(MAKE)
endif

