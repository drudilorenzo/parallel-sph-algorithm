#!/bin/bash

MAX_N=20000

# The one below is the one to use if Hyperthreading is supported
# tot_cores=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')

tot_cores=$(grep -c ^processor /proc/cpuinfo)

echo $tot_cores
for (( core=1; core<=$tot_cores; core++ ))
do
    echo "${core} core" >> results.txt
    # OMP_NUM_THREADS=${core} ./omp-sph 10000 >> results.txt
    # OMP_NUM_THREADS=${core} ./omp-sph 20000 >> results.txt
    for (( n=500; n<=$MAX_N; n+=500 ))
    do
        OMP_NUM_THREADS=${core} ./omp-sph ${n} >> results.txt
    done
done