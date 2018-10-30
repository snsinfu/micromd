DBGFLAGS = \
  -g \
  -fsanitize=address

OPTFLAGS = \
  -Og

CXXFLAGS = \
  -std=c++14 \
  -pedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow \
  $(DBGFLAGS) \
  $(OPTFLAGS)

OBJECTS = \
  main.o \
  test_linear_hash.o \
  test_neighbor_searcher.o \
  test_neighbor_list.o

ARTIFACTS = \
  main \
  *.dSYM \
  $(OBJECTS)


.PHONY: unittest gcov clean
.SUFFIXES: .cc

unittest: main
	./main "$(TAGS)"

clean:
	rm -rf $(ARTIFACTS)

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

main.o: catch.hpp
test_linear_hash.o: catch.hpp ../include/md.hpp
test_neighbor_searcher.o: catch.hpp ../include/md.hpp
test_neighbor_list.o: catch.hpp ../include/md.hpp
