# Parellel programming
This repository contains codes for my assignment programs for the Parallel Programming course at Åbo Akademi University.

## 1. Message-Passing in a Ring using MPI

### Overview
This project focuses on implementing message-passing in a ring using the Message Passing Interface (MPI) in a parallel computing environment. The objective is to facilitate communication between processes arranged in a ring topology, showcasing the power of MPI for parallel programming.

### Implementation
The implementation involves distributing a message in a circular manner, passing it from one process to another until it completes a full loop. This project explores MPI functionalities to demonstrate efficient communication and coordination among distributed computing nodes.

## 2. N-Body Simulation using MPI

### Overview
This project delves into the N-Body simulation problem using the Message Passing Interface (MPI) for parallelization. N-Body simulations model the gravitational interactions between a system of particles, and parallelizing this computation with MPI enhances performance.

### Implementation
The implementation involves dividing the simulation workload across multiple MPI processes, enabling concurrent computation of particle interactions. The project showcases how MPI facilitates parallelism in scientific simulations, improving the efficiency of N-Body simulations.

## 3. Sieve of Eratosthenes using OpenMP

### Overview
The Sieve of Eratosthenes is an algorithm for finding all prime numbers up to a given limit. This project explores parallelization using OpenMP, employing two distinct methods: bit-wise parallelization and byte-wise parallelization.

### Implementation
The project demonstrates how OpenMP directives can be applied to parallelize the Sieve of Eratosthenes algorithm. It includes a bit-wise approach, where individual bits are marked as prime, and a byte-wise approach, where bytes are used for efficient parallelization. The goal is to showcase the benefits of parallel computing in prime number generation.

## How to Run

1. **Message-Passing in a Ring using MPI**
   - Compile the MPI program using the MPI compiler.
   - Run the executable specifying the number of processes to simulate the ring communication.

```bash
cd a1
module load OpenMPI
mpicc ring.c -o ring_mpi
mpirun -np 4 ./ring_mpi
```

2. **N-Body Simulation using MPI**
   - Compile the MPI program for N-Body simulation.
   - Run the executable with the desired parameters, i.e. number of threads

```bash
cd a2
module load OpenMPI

# Sequential N-body sim
mpicc -O2 Nbody.c -o nbody_seq -lm
srun -n 10 ./nbody_seq

# Parallellised N-body sim
mpicc -O2 nbody_par.c -o nbody_par -lm
srun -n 10 ./nbody_par
```

3. **Sieve of Eratosthenes using OpenMP**
   - Compile the OpenMP program for Sieve of Eratosthenes.
   - Run the executable

```bash
cd a3

# Sequential sieve
gcc -fopenmp -O2 seq_sieve.c -o seq_sieve -lm
./seq_sieve

# Parallelized sieve
gcc -fopenmp -O2 par_sieve.c -o par_sieve -lm
./par_sieve
```
