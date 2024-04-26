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

all: epsilon_graph dist_epsilon_graph create_points

install: epsilon_graph
	cp epsilon_graph /global/homes/g/gabeh98/software/cover_tree

.PHONY: version.h

version.h:
	@echo "#define GIT_COMMIT \"$(shell git describe --always --dirty --match 'NOT A TAG')\"" > include/version.h

epsilon_graph: programs/epsilon_graph.cpp src/CoverTree.cpp include/CoverTree.h src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o epsilon_graph $(INCLUDES) $(MPIFLAGS) programs/epsilon_graph.cpp src/CoverTree.cpp src/Point.cpp src/read_args.cpp

dist_epsilon_graph: programs/dist_epsilon_graph.cpp src/DistCoverTree.cpp include/DistCoverTree.h include/Point.h src/MPITimer.cpp include/MPITimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o dist_epsilon_graph $(INCLUDES) $(MPIFLAGS) programs/dist_epsilon_graph.cpp src/DistCoverTree.cpp src/Point.cpp src/MPITimer.cpp src/read_args.cpp

create_points: programs/create_points.cpp src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o create_points $(INCLUDES) $(MPIFLAGS) programs/create_points.cpp src/Point.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs epsilon_graph dist_epsilon_graph create_points version.h
