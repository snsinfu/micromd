OPTFLAGS = \
  -O2 \
  -march=native \
  -mtune=native \
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
  -Wno-c99-extensions \
  $(OPTFLAGS) \
  $(INCLUDES)


.PHONY: all run clean

all: main
	@:

run: main
	./main

clean:
	rm -rf main
