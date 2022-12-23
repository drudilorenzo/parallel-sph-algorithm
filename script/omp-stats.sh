#!/bin/bash

# Threads NOT cores
tot_threads=$(grep -c ^processor /proc/cpuinfo)
input_values=(500 5000 10000)
output_file="omp-stats.txt"
stats_dir="stats"
[ -d ${stats_dir} ] || mkdir -p ${stats_dir}

echo "TOT THREADS: ${tot_threads}" > ${stats_dir}/${output_file}
for (( th=1; th<=${tot_threads}; th++ ))
do
    echo "${th} threads"

    for i in "${input_values[@]}"
    do
        echo "${i}"

        OMP_NUM_THREADS=${th} ./build/sph.omp ${i}
    done
    echo -en "\n"
done >> ${stats_dir}/${output_file}