all: test_point test_brute_force test_cover_tree time_cover_tree

FLAGS=-std=c++17 -O2 -fopenmp

test_point: test_point.cpp Point.o
	g++-13 $(FLAGS) -o $@ $^

test_brute_force: test_brute_force.cpp BruteForce.o Point.o read_args.o
	g++-13 $(FLAGS) -o $@ $^

test_cover_tree: test_cover_tree.cpp CoverTree.o Point.o read_args.o
	g++-13 $(FLAGS) -o $@ $^

time_cover_tree: time_cover_tree.cpp CoverTree.o Point.o read_args.o
	g++-13 $(FLAGS) -o $@ $^

Point.o: Point.cpp Point.h
	g++-13 $(FLAGS) -c -o $@ $<

BruteForce.o: BruteForce.cpp BruteForce.h
	g++-13 $(FLAGS) -c -o $@ $<

CoverTree.o: CoverTree.cpp CoverTree.h
	g++-13 $(FLAGS) -c -o $@ $<

read_args.o: read_args.cpp read_args.h
	g++-13 $(FLAGS) -c -o $@ $<

clean:
	rm -rf a.out *.dSYM *.o
