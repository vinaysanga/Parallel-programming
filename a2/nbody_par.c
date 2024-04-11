// Compile with  gcc -O2 Nbody.c -o Nbody -lm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

const double G = 6.67259e-7; /* Gravitational constant (should be e-10 but modified to get more action */
const double dt = 1.0;       /* Length of timestep */

/* Writes out positions (x,y) of N particles to the file fn
   Returns zero if the file couldn't be opened, otherwise 1 */
int write_particles(int N, double *X, double *Y, char *fn)
{
  FILE *fp;
  /* Open the file */
  if ((fp = fopen(fn, "w")) == NULL)
  {
    printf("Couldn't open file %s\n", fn);
    return 0;
  }
  /* Write the positions to the file fn */
  for (int i = 0; i < N; i++)
  {
    fprintf(fp, "%3.2f %3.2f \n", X[i], Y[i]);
  }
  fprintf(fp, "\n");
  fclose(fp); /* Close the file */
  return (1);
}

// Distance between points with coordinates (px,py) and (qx,qy)
double dist(double px, double py, double qx, double qy)
{
  return sqrt(pow(px - qx, 2) + pow(py - qy, 2));
  // Could also be written as sqrt( (px-qx)*(px-qx) + (py-qy)*(py-qy) )
}

/* Computes forces between bodies */
void ComputeForceParallel(int N, int N_local, double *X, double *Y, double *X_local, double *Y_local, double *mass, double *mass_local, double *Fx, double *Fy, int rank)
{
  const double mindist = 0.0001; /* Minimal distance of two bodies of being in interaction*/

  for (int i = 0; i < N_local; i++)
  {                      // Compute the force for all bodies
    Fx[i] = Fy[i] = 0.0; // Initialize force vector to zero
    for (int j = 0; j < N; j++)
    { // The force on a body i is the sum of forces from all other bodies j
      if ((i + rank * N_local) != j)
      { //     but not from it self
        // Distance between points i and j
        double r = dist(X_local[i], Y_local[i], X[j], Y[j]);

        if (r > mindist)
        {                        // Very near-distance forces are ignored
          double r3 = pow(r, 3); // Could also be written as r3=r*r*r;
          Fx[i] += G * mass_local[i] * mass[j] * (X[j] - X_local[i]) / r3;
          Fy[i] += G * mass_local[i] * mass[j] * (Y[j] - Y_local[i]) / r3;
        }
      }
    }
  }
}

int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);
  int rank, process_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_size);
  const int N = 1000; // Number of bodies
  const int N_local = N / process_size;
  const int timesteps = 1000; // Number of timesteps
  const double size = 100.0;  // Initial positions are in the range [0,100]

  double *mass;
  double *mass_local;
  double *X; /* x-positions of bodies */
  double *Y;
  double *X_local;
  double *Y_local;  /* y-positions of bodies */
  double *Vx_local; /* velocities on x-axis of bodies */
  double *Vy_local; /* velocities on y-axis of bodies */
  double *Fx_local; /* forces on x-axis of bodies */
  double *Fy_local;
  double startTime; /* forces on y-axis of bodies */

  if (rank == 0)
  {
    printf("N-body simulation, number of bodies = %d \n", N);
    fflush(stdout);
  }

  /* Allocate space for variables  */
  mass = (double *)calloc(N, sizeof(double)); // Mass
  X = (double *)calloc(N, sizeof(double));    // Position (x,y) at current time step
  Y = (double *)calloc(N, sizeof(double));
  mass_local = (double *)calloc(N_local, sizeof(double)); // Velocities
  X_local = (double *)calloc(N_local, sizeof(double));    // Velocities
  Y_local = (double *)calloc(N_local, sizeof(double));
  Vx_local = (double *)calloc(N_local, sizeof(double)); // Velocities
  Vy_local = (double *)calloc(N_local, sizeof(double));
  Fx_local = (double *)calloc(N_local, sizeof(double)); // Forces
  Fy_local = (double *)calloc(N_local, sizeof(double));

  // Seed the random number generator so that it generates a fixed sequence
  unsigned short int seedval[3] = {7, 7, 7};
  seed48(seedval);

  /* Initialize mass and position of bodies */

  if (rank == 0)
  {
    for (int i = 0; i < N; i++)
    {
      mass[i] = 1000.0 * drand48(); // 0 <= mass < 1000
      X[i] = size * drand48();      // 0 <= X < 100
      Y[i] = size * drand48();      // 0 <= Y < 100
    }
    // Write intial particle coordinates to a file
    write_particles(N, X, Y, "initial_pos_parallel.txt");
  }

  // send data to local arrays

  MPI_Scatter(X, N_local, MPI_DOUBLE, X_local, N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(Y, N_local, MPI_DOUBLE, Y_local, N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(mass, N_local, MPI_DOUBLE, mass_local, N_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // share whole X Y arrays with all the processes
  MPI_Bcast(X, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(Y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(mass, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    startTime = MPI_Wtime();
  }

  // Compute the initial forces that we get
  ComputeForceParallel(N, N_local, X, Y, X_local, Y_local, mass, mass_local, Fx_local, Fy_local, rank);

  // Set up the velocity vectors caused by initial forces for Leapfrog method
  for (int i = 0; i < N_local; i++)
  {
    Vx_local[i] = 0.5 * dt * Fx_local[i] / mass_local[i];
    Vy_local[i] = 0.5 * dt * Fy_local[i] / mass_local[i];
  }

  /* Main loop:
     - Move the bodies
     - Calculate forces of the bodies with their new position
     - Calculate velocities of the bodies with the new forces
     - Copy the updated positions to the old positions (for use in next timestep)
   */
  int t = 0;
  while (t < timesteps)
  { // Loop for this many timesteps
    t++;
    // Calculate new positions
    for (int i = 0; i < N_local; i++)
    {
      X_local[i] = X_local[i] + Vx_local[i] * dt;
      Y_local[i] = Y_local[i] + Vy_local[i] * dt;
    }

    // Allgather updated positions
    MPI_Allgather(X_local, N_local, MPI_DOUBLE, X, N_local, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(Y_local, N_local, MPI_DOUBLE, Y, N_local, MPI_DOUBLE, MPI_COMM_WORLD);

    /* Calculate forces for the new positions */
    ComputeForceParallel(N, N_local, X, Y, X_local, Y_local, mass, mass_local, Fx_local, Fy_local, rank);

    /* Update velocities of bodies */
    for (int i = 0; i < N_local; i++)
    {
      Vx_local[i] = Vx_local[i] + Fx_local[i] * dt / mass_local[i];
      Vy_local[i] = Vy_local[i] + Fy_local[i] * dt / mass_local[i];
    }

  } /* end of while-loop */

  if (rank == 0)
  {
    printf("\n");
    printf("Time: %6.2f seconds\n", ((MPI_Wtime() - startTime)));
    // Write final particle coordinates to a file
    write_particles(N, X, Y, "final_pos_parallel.txt");
  }
  MPI_Finalize();

  exit(0);
}
