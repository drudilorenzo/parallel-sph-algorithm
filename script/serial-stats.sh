#!/bin/bash

input_values=(500 5000 10000)
output_file="../stats/serial-stats.txt"

for i in "${input_values[@]}"
do
    echo "${i}"

    ../build/sph.serial ${i}

    echo -en "\n"
done > ${output_file}