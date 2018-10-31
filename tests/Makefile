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
  test_forcefield.o \
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
test_linear_hash.o:       catch.hpp ../include/md/detail/neighbor_searcher.hpp
test_neighbor_searcher.o: catch.hpp ../include/md/detail/neighbor_searcher.hpp
test_neighbor_list.o:     catch.hpp ../include/md/detail/neighbor_list.hpp
test_forcefield.o:        catch.hpp ../include/md/forcefield.hpp
test_system.o:            catch.hpp ../include/md/system.hpp
