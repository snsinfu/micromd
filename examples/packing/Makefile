CXX = clang++-brew

OPTFLAGS = \
  -O2 \
  -march=bdver2 \
  -mtune=bdver2 \
  -msse4

INCLUDES = \
  -I ../../include

CXXFLAGS = \
  -std=c++14 \
  -pedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow \
  $(OPTFLAGS) \
  $(INCLUDES)


.PHONY: run clean

run: main
	./main

clean:
	rm -rf main main.dSYM

main:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ main.cc

depends:
	$(CXX) $(CXXFLAGS) -MM -MF- -MT main main.cc > depends.mk

-include depends.mk
