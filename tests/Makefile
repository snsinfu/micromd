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

MDDIR = ../include/md

OBJECTS = \
  main.o \
  test_system.o \
  system/test_attribute_table.o \
  system/detail/test_sum_forcefield.o \
  system/detail/test_dynarray.o \
  forcefield/test_neighbor_pair_forcefield.o \
  forcefield/test_composite_forcefield.o \
  forcefield/detail/test_linear_hash.o \
  forcefield/detail/test_neighbor_searcher.o \
  forcefield/detail/test_neighbor_list.o

ARTIFACTS = \
  main \
  *.dSYM \
  $(OBJECTS)


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

# FIXME: These prerequisites are partial and can cause UB due to ODR violation!

test_system.o: $(MDDIR)/system.hpp

system/test_attribute_table.o :      $(MDDIR)/system/attribute_table.hpp
system/detail/test_sum_forcefield.o: $(MDDIR)/system/detail/sum_forcefield.hpp
system/detail/test_dynarray.o:       $(MDDIR)/system/detail/dynarray.hpp

forcefield/detail/test_linear_hash.o:       $(MDDIR)/forcefield/detail/linear_hash.hpp
forcefield/detail/test_neighbor_searcher.o: $(MDDIR)/forcefield/detail/neighbor_searcher.hpp
forcefield/detail/test_neighbor_list.o:     $(MDDIR)/forcefield/detail/neighbor_list.hpp

forcefield/test_neighbor_pair_forcefield.o: $(MDDIR)/forcefield/neighbor_pair_forcefield.hpp
forcefield/test_composite_forcefield.o:     $(MDDIR)/forcefield/composite_forcefield.hpp
