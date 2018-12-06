DBGFLAGS = \
  -g \
  -fsanitize=address

OPTFLAGS = \
  -Og

INCLUDES = \
  -isystem include \
  -I ../include

override CXXFLAGS += \
  -std=c++14 \
  -pedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow \
  $(DBGFLAGS) \
  $(OPTFLAGS) \
  $(INCLUDES)

SOURCES = $(shell find * -name "*.cc")
OBJECTS = $(SOURCES:.cc=.o)

ARTIFACTS = $(OBJECTS) *.dSYM main


.PHONY: unittest clean depends
.SUFFIXES: .cc

unittest: main
	./main "$(TAGS)"

clean:
	rm -rf $(ARTIFACTS)

depends:
	for s in $(SOURCES); do $(CXX) $(CXXFLAGS) -MM -MF- -MT $${s%.cc}.o $${s}; done > depends.mk

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

-include depends.mk
