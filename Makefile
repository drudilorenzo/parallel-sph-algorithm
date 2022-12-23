# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
#
# Available targets:
#
# - sph: builds the non-GUI version (default)
#
# - sph.gui: builds the GUI version
#
# - all: builds both the GUI and non-GUI versions
#
# - clean: clean up
#
# Last modified on 2022-12-23 by Lorenzo Drudi

SHELL=/bin/bash
EXE:=sph sph.gui sph.omp sph.mpi
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
BUILD_DIR=build
SRC_DIR=src

vpath %.c $(SRC_DIR)

all: create-dir $(EXE)

sph: create-dir sph.serial;

gui: create-dir sph.gui

omp: create-dir sph.omp

mpi: create-dir sph.mpi

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

.PHONY: create-dir clean

create-dir: 
	[ -d $(BUILD_DIR) ] || mkdir -p $(BUILD_DIR) 

clean:
	\rm -f $(BUILD_DIR)/* *.o *~