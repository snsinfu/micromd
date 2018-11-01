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
  detail/test_attribute_table.o \
  detail/test_dynarray.o \
  forcefield/test_composite_forcefield.o \
  forcefield/test_neighbor_pair_forcefield.o \
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

detail/test_linear_hash.o:       ../include/md/detail/linear_hash.hpp
detail/test_neighbor_searcher.o: ../include/md/detail/neighbor_searcher.hpp
detail/test_neighbor_list.o:     ../include/md/detail/neighbor_list.hpp
detail/test_sum_forcefield.o:    ../include/md/detail/sum_forcefield.hpp
detail/test_attribute_table.o:   ../include/md/detail/attribute_table.hpp
detail/test_dynarray.o:          ../include/md/detail/dynarray.hpp

forcefield/test_composite_forcefield.o:     ../include/md/forcefield/composite_forcefield.hpp
forcefield/test_neighbor_pair_forcefield.o: ../include/md/forcefield/neighbor_pair_forcefield.hpp

test_system.o: ../include/md/system.hpp
