D?=0

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
COMPILER=g++
else ifeq ($(UNAME_S),Darwin)
COMPILER=g++-13
endif

ifeq ($(D),1)
FLAGS=-fopenmp -O0 -g -DDEBUG -std=c++17
else
FLAGS=-fopenmp -O2 -std=c++17
endif

INCLUDES=-I./include
PROGS=time_cover_tree

all: $(PROGS)

time_cover_tree: programs/time_cover_tree.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

Point.o: src/Point.cpp include/Point.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

BruteForce.o: src/BruteForce.cpp include/BruteForce.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

CoverTree.o: src/CoverTree.cpp include/CoverTree.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

read_args.o: src/read_args.c include/read_args.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf a.out *.dSYM *.o

distclean: clean
	rm $(PROGS)
