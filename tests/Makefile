DBGFLAGS = \
  -g \
  -fsanitize=address

OPTFLAGS = \
  -Og

INCLUDES = \
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

OBJECTS = \
  main.o \
  detail/test_linear_hash.o \
  detail/test_neighbor_searcher.o \
  detail/test_neighbor_list.o \
  detail/test_sum_forcefield.o \
  test_system.o

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
detail/test_linear_hash.o:       catch.hpp ../include/md/detail/linear_hash.hpp
detail/test_neighbor_searcher.o: catch.hpp ../include/md/detail/neighbor_searcher.hpp
detail/test_neighbor_list.o:     catch.hpp ../include/md/detail/neighbor_list.hpp
detail/test_sum_forcefield.o:    catch.hpp ../include/md/detail/sum_forcefield.hpp
test_system.o:                   catch.hpp ../include/md/system.hpp
