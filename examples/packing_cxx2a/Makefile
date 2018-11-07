OPTFLAGS = \
  -O2 \
  -fno-math-errno \
  -funsafe-math-optimizations

INCLUDES = \
  -I ../../include

CXXFLAGS = \
  -std=c++2a \
  -pedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow \
  $(OPTFLAGS) \
  $(INCLUDES)


.PHONY: all run clean

all: main
	@:

run: main
	./main

clean:
	rm -rf main main.dSYM
