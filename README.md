# parallel-sph-algorithm

Parallel implementation of the [Smoothed Particle Hydrodynamics algorithm](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) developed for the High Performance Computing programming project (Academic Year 2022/2023).

Parallel version of Müller's "Particle-Based Fluid Simulation for Interactive Applications".

[Paper](https://matthias-research.github.io/pages/publications/sca03.pdf)

[Repo](https://github.com/cerrno/mueller-sph)

It offers:

- [x] [Serial implementation](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/sph.c) 
- [x] [Serial implementation with gui](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/sph.c)
- [x] [OpenMP implementation](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/omp-sph.c)
- [x] [MPI implementation](https://github.com/LorenzoDrudi/parallel-sph-algorithm/blob/master/src/mpi-sph.c)

### How to Compile

**Run all the commands in the main directory**.

- ```make all```: compile all the versions 
- ```make sph```: compile serial version 
- ```make gui```: compile gui version
- ```make omp```: compile omp version 
- ```make mpi```: compile mpi version

### How to Run

**Max input size: 20000**.

- *Serial*: ```./build/sph.serial ${INPUT_SIZE}``` (e.g. ```./build/sph.serial 500``` to run with 500 particles)
- *GUI*: ```./build/sph.gui ${INPUT_SIZE}``` (e.g. ```./build/sph.gui 500``` to run with 500 particles)
- *OpenMP*: ```OMP_NUM_THREADS=${NUM_THREADS} ./build/sph.omp ${INPUT_SIZE}``` \
  (e.g. ```OMP_NUM_THREADS=4 ./build/sph.omp 500``` to run with 4 threads and 500 particles)
- *MPI*: ```mpirun -n ${NUM_CORES} ./build/sph.mpi ${INPUT_SIZE}``` \
  (e.g. ```mpirun -n 4 ./build/sph.mpi 500``` to run with 4 cores and 500 particles)
  
  
 
  
  
