# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
# Makefile for the MPI version.
#
# Available targets:
#
# - default: builds the MPI version
#
# - clean: clean up

CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm

all: mpi-sph

mpi-sph: mpi-sph.c
	mpicc ${CFLAGS} $< -o $@ ${LDLIBS}

.PHONY: clean

clean:
	\rm -f mpi-sph *.o *~