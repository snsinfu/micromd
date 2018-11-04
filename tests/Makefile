DBGFLAGS = \
  -g \
  -fsanitize=address

OPTFLAGS = \
  -Og

INCLUDES = \
  -isystem include \
  -I ../include

CXXFLAGS = \
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
OBJECTS = $(addsuffix .o, $(basename $(SOURCES)))

ARTIFACTS = $(OBJECTS) *.dSYM main


.PHONY: unittest clean
.SUFFIXES: .cc

unittest: main
	./main "$(TAGS)"

clean:
	rm -rf $(ARTIFACTS)

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

depends.mk: $(SOURCES)
	for src in $(SOURCES); do $(CXX) $(CXXFLAGS) -MM -MF- -MT $${src%.cc}.o $${src}; done > $@

-include depends.mk
