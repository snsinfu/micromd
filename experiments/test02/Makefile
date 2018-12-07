OPTFLAGS = \
  -O2 \
  -march=native

INCLUDES = \
  -isystem ../../include

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

ARTIFACTS = \
  main \
  samples.dat


.PHONY: all clean

samples.dat: main
	./main > $@

clean:
	rm -f $(ARTIFACTS)
