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

all: dist_benchmark_nng dist_greedy_benchmark

install: dist_benchmark_nng dist_greedy_benchmark
	cp dist_benchmark_nng /global/homes/g/gabeh98/software/cover_tree
	cp dist_greedy_benchmark /global/homes/g/gabeh98/software/cover_tree

dist_benchmark_nng: programs/dist_benchmark_nng.cpp src/CoverTree.cpp include/CoverTree.h src/MPITimer.cpp include/MPITimer.h src/read_args.cpp include/read_args.h
	$(MPICXX) -o dist_benchmark_nng $(INCLUDES) $(MPIFLAGS) programs/dist_benchmark_nng.cpp src/CoverTree.cpp src/MPITimer.cpp src/read_args.cpp

dist_greedy_benchmark: programs/dist_greedy_benchmark.cpp src/MPITimer.cpp include/MPITimer.h src/read_args.cpp include/read_args.h
	$(MPICXX) -o dist_greedy_benchmark $(INCLUDES) $(MPIFLAGS) programs/dist_greedy_benchmark.cpp src/MPITimer.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs dist_benchmark_nng dist_greedy_benchmark
