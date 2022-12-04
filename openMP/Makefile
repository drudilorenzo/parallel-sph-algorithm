# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
# Makefile for the openMP version.
#
# Available targets:
#
# - default: builds the openMP version
#
# - clean: clean up

CFLAGS+=-std=c99 -Wall -Wpedantic -fopenmp
LDLIBS=-lm

all: sph

sph: sph.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

.PHONY: clean

clean:
	\rm -f sph *.o *~