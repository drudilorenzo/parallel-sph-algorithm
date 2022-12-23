#!/bin/bash

input_values=(500 5000 10000)
output_file="serial-stats.txt"
stats_dir="stats"
[ -d ${stats_dir} ] || mkdir -p ${stats_dir}

for i in "${input_values[@]}"
do
    echo "${i}"

    ./build/sph.serial ${i}

    echo -en "\n"
done > ${stats_dir}/${output_file}