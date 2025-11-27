SVARP_VERSION := "1.0"
SVARP_UPDATE  := "Nov 24, 2025"
SVARP_DEBUG   := 0
BUILD_DATE    := "$(shell date)"

# ------------------------------------------------------------------
# Check for zlib1g-dev (needed for HTSlib)
# ------------------------------------------------------------------
ZLIB_DEV := $(shell dpkg -s zlib1g-dev >/dev/null 2>&1 && echo yes || echo no)

ifeq ($(ZLIB_DEV),no)
$(warning zlib1g-dev not found. Please install it with:)
$(warning   sudo apt-get install zlib1g-dev)
endif

# ------------------------------------------------------------------
# Dependency directories
# ------------------------------------------------------------------
DEP_DIR      := dep
HTSLIB_DIR   := $(DEP_DIR)/htslib
WFA_DIR      := $(DEP_DIR)/wfa
WTDBG2_DIR   := $(DEP_DIR)/wtdbg2

HTSLIB_VERSION  := 1.17
HTSLIB_TARBALL  := $(DEP_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_LIB      := $(HTSLIB_DIR)/libhts.a

WFA_VERSION     := 2.3.4
WFA_TARBALL     := $(DEP_DIR)/WFA2-lib-$(WFA_VERSION).tar.gz
WFA_LIB         := $(WFA_DIR)/lib/libwfacpp.a

WTDBG2_BIN      := $(WTDBG2_DIR)/wtdbg2

# ------------------------------------------------------------------
# Compiler / linker flags
# ------------------------------------------------------------------
CXX = g++
CXXFLAGS = -O0 -g -Wall -DDEBUG -std=c++17 \
           -I$(HTSLIB_DIR) -I$(WFA_DIR) \
           -DSVARP_VERSION=\"$(SVARP_VERSION)\" \
           -DBUILD_DATE=\"$(BUILD_DATE)\" \
           -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" \
           -DSVAPR_DEBUG=$(SVARP_DEBUG)


LDFLAGS = $(HTSLIB_LIB) $(WFA_LIB) -lz -lpthread

# ------------------------------------------------------------------
# Project build
# ------------------------------------------------------------------
TARGET_EXEC := svarp
BUILD_DIR   := ./build
SRC_DIRS    := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))


.PHONY: all clean clean-libs libs

all: $(BUILD_DIR)/$(TARGET_EXEC)

# Final link (depends on libs so they exist)
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS) $(HTSLIB_LIB) $(WFA_LIB)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# C++ compile rule
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

clean-libs:
	rm -rf $(DEP_DIR)

# ------------------------------------------------------------------
# Aggregate libs target
# ------------------------------------------------------------------
libs: $(HTSLIB_LIB) $(WFA_LIB) $(WTDBG2_BIN)

# ------------------------------------------------------------------
# HTSLIB 1.22.1
# ------------------------------------------------------------------
$(HTSLIB_LIB): $(HTSLIB_TARBALL)
	mkdir -p $(HTSLIB_DIR)
	tar -xvf $(HTSLIB_TARBALL) -C $(HTSLIB_DIR) --strip-components=1
	cd $(HTSLIB_DIR) && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make

$(HTSLIB_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2 -O $(HTSLIB_TARBALL)

# ------------------------------------------------------------------
# WFA2-lib 2.3.4
# ------------------------------------------------------------------
$(WFA_LIB): $(WFA_TARBALL)
	mkdir -p $(WFA_DIR)
	tar -xzf $(WFA_TARBALL) -C $(WFA_DIR) --strip-components=1
	cd $(WFA_DIR) && make clean all

$(WFA_TARBALL):
	mkdir -p $(DEP_DIR)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v$(WFA_VERSION).tar.gz -O $(WFA_TARBALL)

# ------------------------------------------------------------------
# wtdbg2
# ------------------------------------------------------------------
$(WTDBG2_BIN):
	mkdir -p $(DEP_DIR)
	test -d $(WTDBG2_DIR) || git clone https://github.com/ruanjue/wtdbg2 $(WTDBG2_DIR)
	cd $(WTDBG2_DIR) && make
