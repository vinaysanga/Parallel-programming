//Compile with gcc -O2 -fopenmp par-sieve.c -o par_sieve -lm

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdint.h>
 
#define DEBUG 0        /* Set to 1 if you want a lot of output */
#define NR(i) (2*i+1)  /* The number represented by position i */
#define POS(k) (k/2)   /* The position in the array of number k */

int main (int argc, char *argv[]) {
  uint64_t i , N_pos, N, lastprime=0, nr_primes, k, t;
  char *prime=NULL;
  const char unmarked = (char)0;
  const char marked = (char)1;
  double start, stop;
  
  if (argc < 2) {
    printf("Usage:  %s N\n", argv[0]);
    exit(-1);
  }
  N = strtoull(argv[1], NULL, 10);   // Get the value of N from the command line

  /* To store the marks from 3 to N we need (N-3)/2+1 positions */
  N_pos = (N-3)/2+1;

  start = omp_get_wtime();   /* Start measuring time */

  /* Allocate marks for all odd integers from 3 to N */
  prime = malloc(N_pos);   // Use one byte for each mark
  if (prime == NULL) {
    printf("Could not allocate %d chars of memory\n", N_pos);
    exit(-1);
  }

  /* Mark primes[0] since that is not used */
  prime[0] = marked;
  /* Initialize all odd numbers to unmarked */
  for (i=1; i<N_pos; i++)
    prime[i] = unmarked;

  /* Position i in the array prime now corresponds to the number 2*i+1 */
  for (i=1; NR(i) <= (int)sqrt((double)N); i++) {
    if (prime[i]==unmarked) {   /* Next unmarked position */
      if (DEBUG) printf("Marking multiples of %d: ", NR(i));
      t = NR(i);  /* Position i corresponds to the number t */
      #pragma omp parallel for
      for (k=POS(t*t); k<=N_pos; k+=t) {
        prime[k] = marked;  /* Mark the multiples of i */
        if (DEBUG) printf("%d ", NR(k));
      }
      if (DEBUG) printf("\n");
    }
  }
  
  nr_primes = 1;  /* Remember to count 2 as a prime */
  /* Count the marked numbers */
  #pragma omp parallel default(shared)
  {
    #pragma omp for reduction(+:nr_primes) reduction(max:lastprime)
    for (i=1; i<=N_pos; i++) {
      if (prime[i]==unmarked) {
        lastprime = NR(i);
        nr_primes++;
      }
    }
  }
  
  stop = omp_get_wtime(); 
  printf("Time: %6.2f s\n", (float)(stop-start));

  if (DEBUG) {
    printf("\nPrime numbers smaller than or equal to %d are\n", N);
    printf("%d ", 2);     /* Remember to print the value 2 */
    for (i=1; i<=N_pos; i++) {
      if (prime[i]==unmarked) {
	printf("%d ",NR(i));
      }
    }
  }

  printf("\n%llu primes smaller than or equal to %llu\n", nr_primes, N);
  printf("The largest of these primes is %llu\n", lastprime);

  exit(0);
}
