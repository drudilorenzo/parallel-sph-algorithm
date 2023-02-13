# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
#
#
# Available targets:
#
# - all: builds all the versions (default)
#
# - sph: builds the non-GUI serial version
#
# - gui: builds the GUI serial version
#
# - omp: builds the OpenMP version
#
# - mpi: builds the MPI version
#
# - stats: builds all and run stats scripts (see script folder)
# 
# - clean: clean up (not stats files)
#
# - clean-stats: clean up stats (not executables)
#
#
# Last modified on 2023-01-26 by Drudi Lorenzo - 0000969871

SHELL=/bin/bash
EXE:=omp-sph mpi-sph sph sph.gui
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
STATS_DIR=../stats
SCRIPT_DIR=../script


all: $(EXE)

serial: sph

gui: sph.gui

omp: omp-sph

mpi: mpi-sph

stats: all run-stats 

sph: sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

sph.gui: CFLAGS+=-DGUI
sph.gui: LDLIBS+=-lglut -lGL -lX11
sph.gui: sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

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

