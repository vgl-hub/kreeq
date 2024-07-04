CXX = g++
INCLUDE_DIR = -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = kreeq
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
DECOMPRESSOR_TARGET = decompressor
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LIBS = -lz
LDFLAGS := -pthread

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

SOURCES := main input graph-builder kreeq subgraph kreeq-output variants
OBJECTS := $(addprefix $(BINDIR)/, $(SOURCES))

head: $(OBJECTS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)

debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

all: head validate regenerate decompressor

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h $(GFALIBS_DIR)/include/*.h Makefile | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"

validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp $(LIBS)
	
decompressor: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(DECOMPRESSOR_TARGET) $(SOURCE)/$(DECOMPRESSOR_TARGET).cpp $(LIBS)
	
$(BUILD):
	-mkdir -p $@

$(BINDIR):
	-mkdir -p $@
	
clean:
	$(RM) -r build
	$(MAKE) -C $(GFALIBS_DIR) clean
