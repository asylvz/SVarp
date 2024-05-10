SVARP_VERSION := "1.0"
SVARP_UPDATE := "May 10, 2024"
SVARP_DEBUG := 0
BUILD_DATE := "$(shell date)"

CXX=g++
CXXFLAGS = -O0 -g -Wall -DDEBUG -std=c++17 -I./htslib/ -I./wfa/ -DSVARP_VERSION=\"$(SVARP_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" -DSVAPR_DEBUG=$(SVARP_DEBUG)
LDFLAGS = -lz htslib/libhts.a wfa/lib/libwfacpp.a -lpthread
TARGET_EXEC := svarp
BUILD_DIR := ./build
SRC_DIRS := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)


# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -r $(BUILD_DIR)

libs:
	wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
	mkdir htslib && tar -xvf htslib-1.17.tar.bz2 -C htslib --strip-components=1
	cd htslib && autoconf -i && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make && cd ..
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.4.tar.gz --strip-components=1
	mkdir wfa && tar -xzf v2.3.4.tar.gz -C wfa
	cd wfa && make clean all

.PHONY: clean
