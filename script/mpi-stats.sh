#!/bin/bash

# 
# This script executes the program `mpi-sph` inside the src folder using MPI and
# increasing the cores used from 1 to the number of cores available on your machine.
# Each execution is repeated 5 times.

# The cycle is repeated two times:

# 1. The input size remains constant and the number of cores increases. Execution
# times are useful to compute speedup and strong scaling efficiency.

# 2. The amount of work done by each core remain constant. Execution times are useful 
# to compute weak scaling efficiency.

# Output example: 
# CORES SIZE T1 T2 T3 T4 T5
#
# Change input size and steps to get execution times appropriate to the machine used.
#
# Drudi Lorenzo - 0000969871
# 

# Cores NOT threads
if (( `lscpu | grep Thread | awk '{print $NF}'` == 1  )); then
    CORES=$( grep -c ^processor /proc/cpuinfo )
else
    CORES=$( grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}' )
fi
PROG="../src/mpi-sph"
DEFAULT_SIZE=5000
STEPS=50
REPS=5
STATS_DIR="stats"
STATS_FILE="mpi-stats.txt"

function print_header() {
    echo -en "c\ts\t\t"
    for rep in `seq ${REPS}`; do
        echo -en "t${rep}\t\t"
    done
    echo ""
}

function exec_program() {
    SIZE=$1
    for c in `seq ${CORES}`; do
        echo -en "${c}\t"
        if $2; then
            SIZE=`echo "scale=2; $1 * e(l($c)/2)" | bc -l -q`
        fi
        echo -en "${SIZE}\t"
        for rep in `seq ${REPS}`; do
            EXEC_TIME="$( mpirun -n ${c} ${PROG} ${SIZE} ${STEPS} | tail -1 | sed 's/Elapsed time //' )"
            echo -en "${EXEC_TIME}\t"
        done
        echo ""
    done 
}

if [ ! -f "${PROG}" ]; then
    echo "Program ${PROG} NOT FOUND"
    exit 1
fi

[ -d ../${STATS_DIR} ] || mkdir -p ../${STATS_DIR}

{
    echo "TOT CORES: ${CORES}"

    echo ""
    echo "STRONG SCALING"
    WEAK_SCALING=false
    print_header
    exec_program ${DEFAULT_SIZE} ${WEAK_SCALING}

    echo ""
    echo "WEAK SCALING"
    # Max dim must not exceed 20000. Adapt it on your machine specs
    DEFAULT_SIZE=5000
    WEAK_SCALING=true
    print_header
    exec_program ${DEFAULT_SIZE} ${WEAK_SCALING}
} > ../${STATS_DIR}/${STATS_FILE}