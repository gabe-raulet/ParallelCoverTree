D?=0
DIM?=2
FP?=32
CXX=mpic++
FLAGS=-std=c++20 -I./include -I./misc -DPOINT_DIM=$(DIM) -DFPSIZE=$(FP)

ifeq ($(D),1)
FLAGS+=-O0 -g -fsanitize=address -fno-omit-frame-pointer
else
FLAGS+=-O2
endif

all: create_points build_rgraph graph_diff

create_points: create_points.cpp include/common.h include/ptraits.h include/ptraits.hpp include/timers.h
	$(CXX) -o $@ $(FLAGS) $<

build_rgraph: build_rgraph.cpp include/mpi_comm.h include/mpi_comm.hpp include/bforce.h include/bforce.hpp include/cover_tree.h include/cover_tree.hpp include/common.h include/ptraits.h include/ptraits.hpp include/timers.h
	$(CXX) -o $@ $(FLAGS) $<

graph_diff: graph_diff.cpp include/timers.h
	$(CXX) -o $@ $(FLAGS) $<

clean:
	rm -rf *.dSYM create_points build_rgraph graph_diff

dclean: clean
	git clean -f

