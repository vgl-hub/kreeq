CXX = g++
INCLUDE_DIR = -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = kreeq
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LIBS = -lz
LDFLAGS := -pthread

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

SOURCES := main input kreeq
OBJECTS := $(addprefix $(BINDIR)/, $(SOURCES))

head: $(OBJECTS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)

debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

all: head

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h $(GFALIBS_DIR)/include/*.h Makefile | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"
	
$(BUILD):
	-mkdir -p $@

$(BINDIR):
	-mkdir -p $@
	
clean:
	$(RM) -r build
	$(MAKE) -C $(GFALIBS_DIR) clean
