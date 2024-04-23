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

all: epsilon_graph sgtree

install: epsilon_graph
	cp epsilon_graph /global/homes/g/gabeh98/software/cover_tree

.PHONY: version.h

version.h:
	@echo "#define GIT_COMMIT \"$(shell git describe --always --dirty --match 'NOT A TAG')\"" > include/version.h

epsilon_graph: programs/epsilon_graph.cpp src/CoverTree.cpp include/CoverTree.h src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o epsilon_graph $(INCLUDES) $(MPIFLAGS) programs/epsilon_graph.cpp src/CoverTree.cpp src/Point.cpp src/read_args.cpp

sgtree: programs/sgtree.cpp src/SGTree.cpp include/SGTree.h src/Point.cpp include/Point.h include/MyTimer.h src/read_args.cpp include/read_args.h version.h
	$(MPICXX) -o sgtree $(INCLUDES) $(MPIFLAGS) programs/sgtree.cpp src/SGTree.cpp src/Point.cpp src/read_args.cpp

clean:
	rm -rf *.dSYM *.bin *.fvecs epsilon_graph sgtree version.h
