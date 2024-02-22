FLAGS=-std=c++17 -O2 -fopenmp
INCLUDES=-I./include
PROGS=test_point test_brute_force test_cover_tree time_cover_tree create_dataset

all: $(PROGS)

test_point: programs/test_point.cpp Point.o
	g++-13 $(FLAGS) $(INCLUDES) -o $@ $^

test_brute_force: programs/test_brute_force.cpp BruteForce.o Point.o read_args.o
	g++-13 $(FLAGS) $(INCLUDES) -o $@ $^

test_cover_tree: programs/test_cover_tree.cpp CoverTree.o Point.o read_args.o
	g++-13 $(FLAGS) $(INCLUDES) -o $@ $^

time_cover_tree: programs/time_cover_tree.cpp CoverTree.o Point.o read_args.o
	g++-13 $(FLAGS) $(INCLUDES) -o $@ $^

create_dataset: programs/create_dataset.cpp CoverTree.o Point.o read_args.o
	g++-13 $(FLAGS) $(INCLUDES) -o $@ $^

Point.o: src/Point.cpp include/Point.h
	g++-13 $(FLAGS) $(INCLUDES) -c -o $@ $<

BruteForce.o: src/BruteForce.cpp include/BruteForce.h
	g++-13 $(FLAGS) $(INCLUDES) -c -o $@ $<

CoverTree.o: src/CoverTree.cpp include/CoverTree.h
	g++-13 $(FLAGS) $(INCLUDES) -c -o $@ $<

read_args.o: src/read_args.cpp include/read_args.h
	g++-13 $(FLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf a.out *.dSYM *.o

distclean: clean
	rm $(PROGS)
