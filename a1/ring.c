#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define K 1024            /* One Kilobyte */
#define M K*K             /* One Megabyte */
#define MAXSIZE K*M       /* One Gigabyte */

int main(int argc, char *argv[])
{

    // Initialize MPI
    int np, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Check that we run on more than one process */
    if (np < 2)
    {
        printf("You have to use at least 2 processes to run this program\n");
        MPI_Finalize(); // Quit if there is only one process
        exit(0);
    }

    /* Allocate a sufficiently large message buffer */
    char *buffer = (char *)malloc(MAXSIZE * sizeof(char));

    /* Check that allocation succeeded */
    if (buffer == NULL)
    {
        printf("Could not allocate memory, exiting\n");
        MPI_Finalize();
        exit(0);
    }

    // Initialize the whole message buffer to some values
    for (int i = 0; i < MAXSIZE; i++)
    {
        buffer[i] = (char)75;
    }

    // Repeat until an inputSize <= 0 is entered
    while (1)
    {
        int inputSize;
        if (id == 0)
        {
            printf("Please give an input size in bytes:\n");
            fflush(stdout);
            scanf("%d", &inputSize);
            if (inputSize > MAXSIZE)
            {
                printf("Input size is too large, maximum value is %d\n", MAXSIZE);
                inputSize = 0;
            }
        }

        // Broadcast the message size to all other processes
        MPI_Bcast(&inputSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (inputSize <= 0)
            break; /* Terminate if size = 0 */

        if (id == 0)
        {
            // Send the message to the next process
            double startTime = MPI_Wtime();
            MPI_Send(buffer, inputSize, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(buffer, inputSize, MPI_CHAR, np - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double endTime = MPI_Wtime();
            printf("\nTime taken for communication of message sized %d around the ring with %d processes: %f seconds\n", inputSize, np, endTime - startTime);
        }
        else
        {
            MPI_Recv(buffer, inputSize, MPI_CHAR, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(buffer, inputSize, MPI_CHAR, (id + 1) % np, 0, MPI_COMM_WORLD);
        }
    }

    free(buffer);   // Deallocate the buffer
    MPI_Finalize(); // Terminate MPI
    exit(0);        // Quit
}