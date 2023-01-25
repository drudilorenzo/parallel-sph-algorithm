# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
#
#
# Available targets:
#
# - all: builds all the versions
#
# - omp: builds the OpenMP version
#
# - mpi: builds the MPI version
#
# - stats: builds all and run stats scripts (see script folder)
# 
# - clean: clean up
#
# - clean-stats: clean up stats
#
#
# Last modified on 2023-01-24 by Lorenzo Drudi

SHELL=/bin/bash
EXE:=omp-sph mpi-sph
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
STATS_DIR=../stats
SCRIPT_DIR=../script


all: $(EXE)

omp: omp-sph

mpi: mpi-sph

stats: all run-stats 

omp-sph: CFLAGS+=-fopenmp
omp-sph: omp-sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

mpi-sph: mpi-sph.c
	mpicc ${CFLAGS} $< ${LDLIBS} -o $@

.PHONY: run-stats clean clean-stats

run-stats: 
	$(SCRIPT_DIR)/omp-stats.sh && $(SCRIPT_DIR)/mpi-stats.sh

clean:
	\rm -f $(EXE) *.o *~

clean-stats:
	\rm -f $(STATS_DIR)/*

