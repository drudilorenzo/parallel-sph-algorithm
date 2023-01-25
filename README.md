# parallel-sph-algorithm

Parallel implementation of the [Smoothed Particle Hydrodynamics algorithm](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) developed for the High Performance Computing programming project (Academic Year 2022/2023).

Parallel version of Müller's "Particle-Based Fluid Simulation for Interactive Applications".

[Paper](https://matthias-research.github.io/pages/publications/sca03.pdf)

[Repo](https://github.com/cerrno/mueller-sph)

It offers:

- [x] [OpenMP implementation](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/omp-sph.c)
- [x] [MPI implementation](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/mpi-sph.c)

### How to Compile

**Run all the commands in the main directory**.

- ```make all```: compile all the versions 
- ```make omp```: compile omp version 
- ```make mpi```: compile mpi version
- ```make clean ```: clean up (NOT stats files)

### How to Run

`Max input size`: 20000. \
`Default params`: Particles: 500, Steps: 50.

- *OpenMP*: ```OMP_NUM_THREADS=${NUM_THREADS} omp-sph ${NUM_PARTICLES} ${NUM_STEPS}``` \
  (e.g. ```OMP_NUM_THREADS=4 omp-sph 500``` to run with 4 threads and 500 particles)
- *MPI*: ```mpirun -n ${NUM_CORES} mpi-sph ${NUM_PARTICLES} ${NUM_STEPS}``` \
  (e.g. ```mpirun -n 4 mpi-sph 500``` to run with 4 cores and 500 particles)
  
  
### Stats
  
- ```make stats```: builds all and run the stats scripts (see [script](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/script) folder)
- ```make clean-stats```: clean up stats files

Scripts:
- [omp-stats](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/script/omp-stats.sh)
- [mpi-stats](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/script/mpi-stats.sh) 

Each script executes the program several times increasing the number of threads/cores from 1 to the number available on your machine.
Each execution is repeated 5 times.
The cycle is done two times:
1. The input size remains constant and the number of threads/cores increases. Execution
   times are useful to get speedup and strong scaling efficiency;
2. The amount of work done by each thread/core remain constant. The execution times are useful to get weak scaling efficiency.

