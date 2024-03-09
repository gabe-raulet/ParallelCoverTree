D?=0
UNAME_S:=$(shell uname -s)

ifeq ($(D),1)
FLAGS=-fopenmp -O0 -g -std=c++17 -pedantic
else
FLAGS=-fopenmp -O2 -DNDEBUG -std=c++17
endif

ifeq ($(UNAME_S), Linux)
COMPILER=g++
FLAGS+=-fsanitize=address -fno-omit-frame-pointer
else
COMPILER=g++-13
endif

INCLUDES=-I./include

all: cluster_benchmark build_benchmark create_data

test: cluster_benchmark
	@./cluster_benchmark -i testdata/points.1K.2d.fvecs -r 0.075 -o testdata/points.1K.2d.r075.clusters.compare
	@diff -s testdata/points.1K.2d.r075.clusters testdata/points.1K.2d.r075.clusters.compare
	@./cluster_benchmark -i testdata/points.1K.4d.fvecs -r 0.3 -o testdata/points.1K.4d.r3.clusters.compare
	@diff -s testdata/points.1K.4d.r3.clusters testdata/points.1K.4d.r3.clusters.compare
	@./cluster_benchmark -i testdata/points.1K.8d.fvecs -r 0.975 -o testdata/points.1K.8d.r975.clusters.compare
	@diff -s testdata/points.1K.8d.r975.clusters testdata/points.1K.8d.r975.clusters.compare
	@rm -rf testdata/*.clusters.compare

cluster_benchmark: programs/cluster_benchmark.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(COMPILER) -o cluster_benchmark $(INCLUDES) $(FLAGS) programs/cluster_benchmark.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

build_benchmark: programs/build_benchmark.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(COMPILER) -o build_benchmark $(INCLUDES) $(FLAGS) programs/build_benchmark.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

create_data: programs/create_data.cpp src/read_args.cpp include/read_args.h
	$(COMPILER) -o create_data $(INCLUDES) $(FLAGS) programs/create_data.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs cluster_benchmark build_benchmark create_data
