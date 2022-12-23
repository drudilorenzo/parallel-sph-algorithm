# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
#
# The source files are in the src folder and the executables in the build one.
#
# Available targets:
#
# - all: builds all the versions
#
# - sph: builds the non-GUI version
#
# - gui: builds the GUI version
#
# - omp: builds the openMP version
#
# - mpi: builds the MPI version
#
# - create-build-dir: create the build folder
#
# - stats: run stats scripts (see script folder)
# 
# - clean: clean up
#
# - clean-stats: clean up stats
#
#
# Last modified on 2022-12-23 by Lorenzo Drudi

SHELL=/bin/bash
EXE:=sph.serial sph.gui sph.omp sph.mpi
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
BUILD_DIR=build
SRC_DIR=src
STATS_DIR=stats
SCRIPT_DIR=script

vpath %.c $(SRC_DIR)

all: create-build-dir $(EXE)

sph: create-build-dir sph.serial;

gui: create-build-dir sph.gui

omp: create-build-dir sph.omp

mpi: create-build-dir sph.mpi

stats: clean-stats all run-stats 

sph.serial:: sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $(BUILD_DIR)/$@

sph.gui: CFLAGS+=-DGUI
sph.gui: LDLIBS+=-lglut -lGL -lX11
sph.gui: sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $(BUILD_DIR)/$@

sph.omp: CFLAGS+=-fopenmp
sph.omp: omp-sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $(BUILD_DIR)/$@

sph.mpi: mpi-sph.c
	mpicc ${CFLAGS} $< ${LDLIBS} -o $(BUILD_DIR)/$@

.PHONY: create-dir clean clean-stats run-stats

create-build-dir: 
	[ -d $(BUILD_DIR) ] || mkdir -p $(BUILD_DIR) 

run-stats: 
	./$(SCRIPT_DIR)/serial-stats.sh && ./$(SCRIPT_DIR)/omp-stats.sh && ./$(SCRIPT_DIR)/mpi-stats.sh

clean:
	\rm -f $(BUILD_DIR)/* *.o *~

clean-stats:
	\rm -f $(STATS_DIR)/*

