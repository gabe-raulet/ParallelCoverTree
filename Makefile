D?=0
UNAME_S:=$(shell uname -s)

ifeq ($(D),1)
FLAGS=-fopenmp -O0 -g -std=c++17 -pedantic
MPIFLAGS=-O0 -g -std=c++17 -fsanitize=address -fno-omit-frame-pointer
else
FLAGS=-fopenmp -O2 -std=c++17
MPIFLAGS=-O2 -std=c++17
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

all: create_points build_graph dist_build_graph

install: create_points build_graph dist_build_graph
	cp create_points /global/homes/g/gabeh98/software/cover_tree
	cp build_graph /global/homes/g/gabeh98/software/cover_tree
	cp dist_build_graph /global/homes/g/gabeh98/software/cover_tree

.PHONY: version.h

version.h:
	@echo "#define GIT_COMMIT \"$(shell git describe --always --dirty --match 'NOT A TAG')\"" > include/version.h

create_points: programs/create_points.cpp src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o create_points $(INCLUDES) $(MPIFLAGS) programs/create_points.cpp src/Point.cpp src/read_args.cpp

build_graph: programs/build_graph.cpp src/CoverTree.cpp include/CoverTree.h src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o build_graph $(INCLUDES) $(MPIFLAGS) programs/build_graph.cpp src/CoverTree.cpp src/Point.cpp src/read_args.cpp

dist_build_graph: programs/dist_build_graph.cpp src/DistCoverTree.cpp include/DistCoverTree.h src/CoverTree.cpp include/CoverTree.h src/Point.cpp include/Point.h src/MPITimer.cpp include/MPITimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o dist_build_graph $(INCLUDES) $(MPIFLAGS) programs/dist_build_graph.cpp src/DistCoverTree.cpp src/CoverTree.cpp src/Point.cpp src/MPITimer.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs create_points build_graph dist_build_graph version.h
