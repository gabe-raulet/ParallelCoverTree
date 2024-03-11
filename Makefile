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
MPICXX=g++
FLAGS+=-fsanitize=address -fno-omit-frame-pointer
MPIFLAGS+=-fsanitize=address -fno-omit-frame-pointer
else
CXX=g++-13
MPICXX=mpic++
endif

INCLUDES=-I./include

all: build_tree build_epsilon_graph dist_build_epsilon_graph create_data

install: cluster_test build_tree build_epsilon_graph create_data
	cp cluster_test /global/homes/g/gabeh98/software/cover_tree
	cp build_tree /global/homes/g/gabeh98/software/cover_tree
	cp build_epsilon_graph /global/homes/g/gabeh98/software/cover_tree
	cp create_data /global/homes/g/gabeh98/software/cover_tree

test: cluster_test build_tree
	@./build_tree -i testdata/points.1K.2d.fvecs -o testdata/1K.2d.cover_tree
	@./cluster_test -i testdata/1K.2d.cover_tree -r 0.075 -o testdata/points.1K.2d.r075.clusters.compare
	@diff -s testdata/points.1K.2d.r075.clusters testdata/points.1K.2d.r075.clusters.compare
	@./build_tree -i testdata/points.1K.4d.fvecs -o testdata/1K.4d.cover_tree
	@./cluster_test -i testdata/1K.4d.cover_tree -r 0.3 -o testdata/points.1K.4d.r3.clusters.compare
	@diff -s testdata/points.1K.4d.r3.clusters testdata/points.1K.4d.r3.clusters.compare
	@./build_tree -i testdata/points.1K.8d.fvecs -o testdata/1K.8d.cover_tree
	@./cluster_test -i testdata/1K.8d.cover_tree -r 0.975 -o testdata/points.1K.8d.r975.clusters.compare
	@diff -s testdata/points.1K.8d.r975.clusters testdata/points.1K.8d.r975.clusters.compare
	@rm -rf testdata/*.clusters.compare testdata/*.cover_tree

cluster_test: programs/cluster_test.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(CXX) -o cluster_test $(INCLUDES) $(FLAGS) programs/cluster_test.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

build_tree: programs/build_tree.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(CXX) -o build_tree $(INCLUDES) $(FLAGS) programs/build_tree.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

build_epsilon_graph: programs/build_epsilon_graph.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(CXX) -o build_epsilon_graph $(INCLUDES) $(FLAGS) programs/build_epsilon_graph.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

dist_build_epsilon_graph: programs/dist_build_epsilon_graph.cpp src/CoverTree.cpp include/CoverTree.h src/VectorIO.cpp include/VectorIO.h src/read_args.cpp include/read_args.h
	$(MPICXX) -o dist_build_epsilon_graph $(INCLUDES) $(MPIFLAGS) programs/dist_build_epsilon_graph.cpp src/CoverTree.cpp src/VectorIO.cpp src/read_args.cpp

create_data: programs/create_data.cpp src/read_args.cpp include/read_args.h
	$(CXX) -o create_data $(INCLUDES) $(FLAGS) programs/create_data.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs cluster_test build_tree build_epsilon_graph dist_build_epsilon_graph create_data
