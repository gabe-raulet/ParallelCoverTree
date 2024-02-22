all: test_point test_brute_force test_cover_tree

FLAGS=-std=c++17 -g -O0 -fsanitize=address -fno-omit-frame-pointer

test_point: test_point.cpp Point.o
	clang++ $(FLAGS) -o $@ $^

test_brute_force: test_brute_force.cpp BruteForce.o Point.o read_args.o
	clang++ $(FLAGS) -o $@ $^

test_cover_tree: test_cover_tree.cpp CoverTree.o Point.o read_args.o
	clang++ $(FLAGS) -o $@ $^

Point.o: Point.cpp Point.h
	clang++ $(FLAGS) -c -o $@ $<

BruteForce.o: BruteForce.cpp BruteForce.h
	clang++ $(FLAGS) -c -o $@ $<

CoverTree.o: CoverTree.cpp CoverTree.h
	clang++ $(FLAGS) -c -o $@ $<

read_args.o: read_args.cpp read_args.h
	clang++ $(FLAGS) -c -o $@ $<

clean:
	rm -rf a.out *.dSYM *.o
