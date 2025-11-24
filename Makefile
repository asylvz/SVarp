SVARP_VERSION := "1.0"
SVARP_UPDATE := "Nov 24, 2025"
SVARP_DEBUG := 0
BUILD_DATE := "$(shell date)"

CXX = g++
CXXFLAGS = -O0 -g -Wall -DDEBUG -std=c++17 \
           -I./htslib/ -I./wfa/ \
           -DSVARP_VERSION=\"$(SVARP_VERSION)\" \
           -DBUILD_DATE=\"$(BUILD_DATE)\" \
           -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" \
           -DSVAPR_DEBUG=$(SVARP_DEBUG)

LDFLAGS = -lz htslib/libhts.a wfa/lib/libwfacpp.a -lpthread

TARGET_EXEC := svarp
BUILD_DIR := ./build
SRC_DIRS := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

.PHONY: all clean libs

all: $(BUILD_DIR)/$(TARGET_EXEC)

# Final link
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# C++ compile rule
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

libs:
	# HTSLIB
	wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2
	mkdir -p htslib
	tar -xvf htslib-1.22.1.tar.bz2 -C htslib --strip-components=1
	cd htslib && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make

	# WFA2-lib
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.4.tar.gz
	mkdir -p wfa
	tar -xzf v2.3.4.tar.gz --strip-components=1 -C wfa
	cd wfa && make clean all

