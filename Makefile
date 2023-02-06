SVARP_VERSION := "0.1"
SVARP_UPDATE := "February 1, 2023"
SVARP_DEBUG := 0
BUILD_DATE := "$(shell date)"

CXX=g++
CXXFLAGS = -g -O0 -Wall -DDEBUG -std=c++17 -DSVARP_VERSION=\"$(SVARP_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DSVARP_UPDATE=\"$(SVARP_UPDATE)\" -DSVAPR_DEBUG=$(SVARP_DEBUG)
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

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
