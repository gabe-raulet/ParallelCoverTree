D?=0

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
COMPILER=g++
else ifeq ($(UNAME_S),Darwin)
COMPILER=g++-13
endif

ifeq ($(D),1)
FLAGS=-fopenmp -O0 -g -DDEBUG -std=c++17 -pedantic
else
FLAGS=-fopenmp -O2 -std=c++17
endif

INCLUDES=-I./include
PROGS=test_cover_tree test_queries test_all test_construction create_data

all: $(PROGS)

test_all: programs/test_all.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

test_construction: programs/test_construction.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

test_cover_tree: programs/test_cover_tree.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

create_data: programs/create_data.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

test_queries: programs/test_queries.cpp CoverTree.o Point.o read_args.o
	$(COMPILER) $(FLAGS) $(INCLUDES) -o $@ $^

Point.o: src/Point.cpp include/Point.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

CoverTree.o: src/CoverTree.cpp include/CoverTree.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

read_args.o: src/read_args.c include/read_args.h
	$(COMPILER) $(FLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf a.out *.dSYM *.o *.bin

distclean: clean
	rm -f $(PROGS)
