#!/bin/bash

# Cores NOT threads
tot_cores=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
input_values=(500 5000 10000)
output_file="../stats/mpi-stats.txt"

echo "TOT CORES: ${tot_cores}" > ${output_file}
for (( c=1; c<=${tot_cores}; c++ ))
do
    echo "${c} cores" >> ${output_file}

    for i in "${input_values[@]}"
    do
        echo "${i}"

        mpirun -n ${c} ../build/sph.mpi ${i} 
    done

    echo -en "\n"
done >> ${output_file}