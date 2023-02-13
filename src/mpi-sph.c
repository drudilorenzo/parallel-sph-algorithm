/****************************************************************************
 *
 * sph.c -- Smoothed Particle Hydrodynamics
 *
 * https://github.com/cerrno/mueller-sph
 *
 * Copyright (C) 2016 Lucas V. Schuermann
 * Copyright (C) 2022 Moreno Marzolla
 * Copyright (C) 2023 Lorenzo Drudi
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 ****************************************************************************/

/*
    Student: Drudi Lorenzo - 0000969871

    Compile with:
    mpicc -std=c99 -Wall -Wpedantic mpi-sph.c -o mpi-sph -lm

    Run with:
    mpirun -n ${N_THREADS} ./mpi-sph ${N_PARTICLES} ${N_STEPS}
    Default params: Particles: 500, Steps: 50

    mpirun -n 4 ./mpi-sph
*/

/* It must be the programme's first include.
   It's used to take the wall-clock time. */
#include "hpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stddef.h>
#include <mpi.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* "Particle-Based Fluid Simulation for Interactive Applications" by
   Müller et al. solver parameters */

const float Gx = 0.0, Gy = -10.0;   // external (gravitational) forces
const float REST_DENS = 300;    // rest density
const float GAS_CONST = 2000;   // const for equation of state
const float H = 16;             // kernel radius
const float EPS = 16;           // equal to H
const float MASS = 2.5;         // assume all particles have the same mass
const float VISC = 200;         // viscosity constant
const float DT = 0.0007;        // integration timestep
const float BOUND_DAMPING = -0.5;

const int MAX_PARTICLES = 20000;
// Larger window size to accommodate more particles
#define WINDOW_WIDTH 3000
#define WINDOW_HEIGHT 2000

const int DAM_PARTICLES = 500;

const float VIEW_WIDTH = 1.5 * WINDOW_WIDTH;
const float VIEW_HEIGHT = 1.5 * WINDOW_HEIGHT;

/* Particle data structure; stores position, velocity, and force for
   integration stores density (rho) and pressure values for SPH.

   You may choose a different layout of the particles[] data structure
   to suit your needs. */
typedef struct {
    float x, y;         // position
    float vx, vy;       // velocity
    float fx, fy;       // force
    float rho, p;       // density, pressure
} particle_t;


/* Main particles array.

   It's managed by the master process (rank == 0).
   Each process send its local updated cells to this array. */
particle_t *particles; 

/* Local particles array.

   Each process do the computation on its domain partition and then
   send the local updated particles to the main array. */
particle_t *local_particles; 

int n_particles = 0;    // number of currently active particles
MPI_Datatype mpi_particles; // MPI data type to send/receive data of type particle_t

/* Counts array.

   Entry i specifies the number of particles computed by process i. */
int *counts;

/* Displacements array.

   Entry i specifies the displacement of process i data from the starting 
   of the array. */
int *displs; 

int my_rank; // Rank of the process

/**
 * Return a random value in [a, b]
 */
float randab(float a, float b)
{
    return a + (b-a)*rand() / (float)(RAND_MAX);
}

/**
 * Set initial position of particle `*p` to (x, y); initialize all
 * other attributes to default values (zeros).
 */
void init_particle(particle_t *p, float x, float y)
{
    p->x = x;
    p->y = y;
    p->vx = p->vy = 0.0;
    p->fx = p->fy = 0.0;
    p->rho = 0.0;
    p->p = 0.0;
}

/**
 * Return nonzero if (x, y) is within the frame
 */
int is_in_domain(float x, float y)
{
    return ((x < VIEW_WIDTH - EPS) &&
            (x > EPS) &&
            (y < VIEW_HEIGHT - EPS) &&
            (y > EPS));
}

/**
 * Initialize the SPH model with `n` particles. The caller is
 * responsible for allocating the `particles[]` array of size
 * `MAX_PARTICLES`.
 *
 * DO NOT parallelize this function, since it calls rand() which is
 * not thread-safe.
 */
void init_sph(int n)
{
    n_particles = 0;
    printf("Initializing with %d particles\n", n);

    for (float y = EPS; y < VIEW_HEIGHT - EPS; y += H) {
        for (float x = EPS; x <= VIEW_WIDTH * 0.8f; x += H) {
            if (n_particles < n) {
                float jitter = rand() / (float)RAND_MAX;
                init_particle(particles + n_particles, x+jitter, y);
                n_particles++;
            } else {
                return;
            }
        }
    }
    assert(n_particles == n);
}

void compute_density_pressure(void)
{
    const float HSQ = H * H;    // radius^2 for optimization

    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. 
       
       It modifies only its local particles.*/
    const float POLY6 = 4.0 / (M_PI * pow(H, 8));
    const int local_n = counts[my_rank];
    for (int i = 0; i < local_n; i++) {
        particle_t *pi = &local_particles[i];
        pi->rho = 0.0;
        for (int j = 0; j < n_particles; j++) {
            const particle_t *pj = &particles[j];

            const float dx = pj->x - pi->x;
            const float dy = pj->y - pi->y;
            const float d2 = dx*dx + dy*dy;

            if (d2 < HSQ) {
                pi->rho += MASS * POLY6 * pow(HSQ - d2, 3.0);
            }
        }
        pi->p = GAS_CONST * (pi->rho - REST_DENS);
    }
}

void compute_forces(void)
{
    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. 
       
       Each process modifies only the particles of its partition.*/
    const float SPIKY_GRAD = -10.0 / (M_PI * pow(H, 5));
    const float VISC_LAP = 40.0 / (M_PI * pow(H, 5));
    const float EPS = 1e-6;
    const int local_n = counts[my_rank];

    for (int i = 0; i < local_n; i++) {
        particle_t *pi = &local_particles[i];
        float fpress_x = 0.0, fpress_y = 0.0;
        float fvisc_x = 0.0, fvisc_y = 0.0;

        for (int j = 0; j < n_particles; j++) {
            const particle_t *pj = &particles[j];

            if (pi == pj)
                continue;

            const float dx = pj->x - pi->x;
            const float dy = pj->y - pi->y;
            const float dist = hypotf(dx, dy) + EPS; // avoids division by zero later on

            if (dist < H) {
                const float norm_dx = dx / dist;
                const float norm_dy = dy / dist;
                // compute pressure force contribution
                fpress_x += -norm_dx * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
                fpress_y += -norm_dy * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
                // compute viscosity force contribution
                fvisc_x += VISC * MASS * (pj->vx - pi->vx) / pj->rho * VISC_LAP * (H - dist);
                fvisc_y += VISC * MASS * (pj->vy - pi->vy) / pj->rho * VISC_LAP * (H - dist);
            }
        }
        const float fgrav_x = Gx * MASS / pi->rho;
        const float fgrav_y = Gy * MASS / pi->rho;
        pi->fx = fpress_x + fvisc_x + fgrav_x;
        pi->fy = fpress_y + fvisc_y + fgrav_y;
    }
}

void integrate(void)
{
    /* Each process modifies only the particles of its partition. */
    const int local_n = counts[my_rank];
    for (int i = 0; i < local_n; i++) {
        particle_t *p = &local_particles[i];
        // forward Euler integration
        p->vx += DT * p->fx / p->rho;
        p->vy += DT * p->fy / p->rho;
        p->x += DT * p->vx;
        p->y += DT * p->vy;

        // enforce boundary conditions
        if (p->x - EPS < 0.0) {
            p->vx *= BOUND_DAMPING;
            p->x = EPS;
        }
        if (p->x + EPS > VIEW_WIDTH) {
            p->vx *= BOUND_DAMPING;
            p->x = VIEW_WIDTH - EPS;
        }
        if (p->y - EPS < 0.0) {
            p->vy *= BOUND_DAMPING;
            p->y = EPS;
        }
        if (p->y + EPS > VIEW_HEIGHT) {
            p->vy *= BOUND_DAMPING;
            p->y = VIEW_HEIGHT - EPS;
        }
    }
}

float avg_velocities( void )
{
    double local_result = 0.0;
    const int local_n = counts[my_rank];

    for (int i = 0; i < local_n; i++) {
        /* the hypot(x,y) function is equivalent to sqrt(x*x + y*y); */
        local_result += hypot(local_particles[i].vx, local_particles[i].vy) / n_particles;
    }

    return local_result;
}

void update( void )
{
    compute_density_pressure();

    /* After the `compute_density_pressure` function we have to gather 
       the updated values and then distribute them to all the processes
       since are needed by the next step. */
    MPI_Allgatherv(
        local_particles,
        counts[my_rank],
        mpi_particles,
        particles,
        counts,
        displs,
        mpi_particles,
        MPI_COMM_WORLD
    );

    compute_forces();
    integrate();

    /* As we've done before we gather and distribute all the updated particles. */
    MPI_Allgatherv(
        local_particles,
        counts[my_rank],
        mpi_particles,
        particles,
        counts,
        displs,
        mpi_particles,
        MPI_COMM_WORLD
    );
}

int main(int argc, char **argv)
{
    int n = DAM_PARTICLES;
    int nsteps = 50;
    int comm_sz;

    srand(1234);
 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    if (argc > 3) {
        fprintf(stderr, "Usage: %s [nparticles [nsteps]]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (argc > 1) {
        n = atoi(argv[1]);
    }

    if (argc > 2) {
        nsteps = atoi(argv[2]);
    }

    if (n > MAX_PARTICLES) {
        fprintf(stderr, "FATAL: the maximum number of particles is %d\n", MAX_PARTICLES);
        return EXIT_FAILURE;
    }

    particles = (particle_t*)malloc(n * sizeof(*particles));
    assert( particles != NULL );
    if (my_rank == 0) {
        init_sph(n);
    }

    /* Compute the counts and displacements array.
       They're used to partition domains of arbitrary size. */
    counts = (int*)malloc(comm_sz * sizeof(*counts)); assert(counts != NULL);
    displs = (int*)malloc(comm_sz * sizeof(*displs)); assert(displs != NULL);
    for (int i = 0; i < comm_sz; i++) {
        const int start = n * i / comm_sz;
        const int end = n * (i+1) / comm_sz;
        counts[i] = end - start;
        displs[i] = start;
    }

    /* Local array used to store the partition of each process. */
    const int local_n = counts[my_rank];
    local_particles = (particle_t*)malloc(local_n * sizeof(*local_particles)); 
    assert(local_particles != NULL);

    /* MPI type struct to send/receive data of type particles_t.
       MPI already permits to send/receive contiguous elements without 
       a custom data type. Despite that I prefered to use it. */
    MPI_Type_contiguous(8, MPI_FLOAT, &mpi_particles);
    MPI_Type_commit(&mpi_particles);

    /* Broadcast the particles array to each processor. */
    MPI_Bcast(
        particles,
        n,
        mpi_particles,
        0,
        MPI_COMM_WORLD
    ); 

    /* Since particles is a dynamic array, it's not possible to get its size.
       So we need to broadcast also its length. */
    MPI_Bcast(
        &n_particles,
        1,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    /* Since we've already sent all the particles array, is not needed a scatterv
       to partition the domain into the local arrays.
       We can do that just with a for cycle using the process' displacement. */
    for (int i = 0; i < local_n; i++) {
        const int global_index = displs[my_rank] + i; 
        local_particles[i] = particles[global_index];
    }

    const float tstart = hpc_gettime();
    for (int s=0; s<nsteps; s++) {
        update();

        /* the average  velocities MUST be computed at each step, even
           if it is not shown (to ensure constant workload per
           iteration) */
        const float local_avg = avg_velocities();
        float avg = 0;
        MPI_Allreduce(
            &local_avg,
            &avg,
            1,
            MPI_FLOAT,
            MPI_SUM,
            MPI_COMM_WORLD
        );

        if (my_rank == 0 && s % 10 == 0) {
            printf("step %5d, avgV=%f\n", s, avg);
        }
    }
    const float elapsed = hpc_gettime() - tstart;
    if (my_rank == 0) {        
        printf("Elapsed time %.2f\n", elapsed);
    }

    free(particles);
    free(local_particles);
    free(counts);
    free(displs);

    MPI_Type_free(&mpi_particles);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
