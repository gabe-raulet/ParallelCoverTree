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

all: benchmark_nng benchmark_nng_bf dist_benchmark_nng_bf

install: benchmark_nng benchmark_nng_bf dist_benchmark_nng_bf
	cp benchmark_nng /global/homes/g/gabeh98/software/cover_tree
	cp benchmark_nng_bf /global/homes/g/gabeh98/software/cover_tree
	cp dist_benchmark_nng_bf /global/homes/g/gabeh98/software/cover_tree

benchmark_nng: programs/benchmark_nng.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(CXX) -o benchmark_nng $(INCLUDES) $(FLAGS) programs/benchmark_nng.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

benchmark_nng_bf: programs/benchmark_nng_bf.cpp src/read_args.cpp include/read_args.h
	$(CXX) -o benchmark_nng_bf $(INCLUDES) $(FLAGS) programs/benchmark_nng_bf.cpp src/read_args.cpp

dist_benchmark_nng_bf: programs/dist_benchmark_nng_bf.cpp src/read_args.cpp include/read_args.h
	$(MPICXX) -o dist_benchmark_nng_bf $(INCLUDES) $(MPIFLAGS) programs/dist_benchmark_nng_bf.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs benchmark_nng benchmark_nng_bf dist_benchmark_nng_bf
