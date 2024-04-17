D?=0
UNAME_S:=$(shell uname -s)

ifeq ($(D),1)
FLAGS=-fopenmp -O0 -g -std=c++17 -pedantic
MPIFLAGS=-O0 -g -std=c++17
else
FLAGS=-fopenmp -O2 -DNDEBUG -std=c++17
MPIFLAGS=-O2 -DNDEBUG -std=c++17
endif

ifeq ($(UNAME_S), Linux)
CXX=g++
MPICXX=mpic++
FLAGS+=-fsanitize=address -fno-omit-frame-pointer
else
CXX=g++-13
MPICXX=mpic++
endif

INCLUDES=-I./include

all: epsilon_graph

install: epsilon_graph
	cp epsilon_graph /global/homes/g/gabeh98/software/cover_tree

epsilon_graph: programs/epsilon_graph.cpp src/CoverTree.cpp include/CoverTree.h src/Point.cpp include/Point.h src/read_args.cpp include/read_args.h
	$(CXX) -o epsilon_graph $(INCLUDES) $(MPIFLAGS) programs/epsilon_graph.cpp src/CoverTree.cpp src/Point.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs epsilon_graph
