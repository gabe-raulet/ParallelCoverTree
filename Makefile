all: test_point

FLAGS=-std=c++17 -g -O0 -fsanitize=address -fno-omit-frame-pointer

test_point: test_point.cpp Point.o
	clang++ $(FLAGS) -o test_point test_point.cpp Point.o

Point.o: Point.cpp Point.h
	clang++ $(FLAGS) -c -o Point.o Point.cpp

clean:
	rm -rf a.out *.dSYM *.o
