#include "proj_common.h"
#include <stdio.h>
#include <omp.h>
#include <string.h>

/**
 * @brief Main function of the simulation.
 *
 * This function performs the entire simulation loop for the Vicsek model of Flocking Birds using MPI Processes and OpenMP Threads.
 * Using functions in proj_common.h the model is solved and results are printed in one line per timestep. 
 * The resulting values can then be pasted into the Python program to visualize the flock.
 *
 * @param argc Number of command-line arguments. Unused.
 * @param argv Array of command-line argument strings. Unused.
 * @return Returns 0 upon successful completion.
 */
int main(int argc, char *argv[])
{
    int i, j, rank, size, provided, num_pp, startnum; /**< Loop counters and MPI variables. */
    struct Bird *b; /**< Reusable pointer to a bird struct. */

    remove("output.dat");

    double R = pow(R_INIT, 2); /**< Interaction radius squared, used in Pythagorean theorem. */

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided); /**< Initialize MPI with single thread support. */
    MPI_Comm_size(MPI_COMM_WORLD, &size); /**< Get the total number of processes. */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /**< Get the rank of the current process. */
    
    MPI_File fh; /**< Creates a MPI IO File Reference*/
    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); /**< Opens the MPI IO File as write only*/

    if (rank == 0) {
        if (NUMBER % size != 0) { //**< Program will only run if amount of processes is a divisor of total birds for maths reasons. */
            printf("Please start with a process number that is a divisor of %d", NUMBER);
            return 1;
        }
        printf("%d %d %f %f\n", NUMBER, TIMESTEPS, L, DT);
    } 

    num_pp = NUMBER / size; /**< Calculate the number of birds per process. */
    startnum = rank * num_pp; /**< Calculate the starting index for birds for this process. */
    
    srand(SEED * rank); /**< Seed the random number generator with a unique seed for each process. */

    struct Bird *birds = calloc(NUMBER, sizeof(struct Bird)); /**< Allocate memory for all birds. */
    struct Bird *proc_birds = calloc(num_pp, sizeof(struct Bird)); /**< Allocate memory for birds of this process. */
    double *print_values = malloc(num_pp * PRINT_DOUBLES * sizeof(double)); /**< Allocate memory for values to save from this process. */

    double startTime = omp_get_wtime(); /**< Record start time of simulation. */

    if (rank == 0) {
        #pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < NUMBER; i++) { /**< Initialize birds only on process 0. */
            initBird(&birds[i]);
        }
    }
    MPI_Bcast(birds, NUMBER * BIRD_DOUBLES, MPI_DOUBLE, 0, MPI_COMM_WORLD); /**< Broadcast all bird data to all processes. */

    memcpy(proc_birds, &birds[startnum], sizeof(struct Bird) * num_pp); /**< Copy the birds for this process. */

    for (i = 0; i < TIMESTEPS; i++) { /**< Main simulation loop. */
        #pragma omp parallel for schedule(static) private(j)
        for (j = 0; j < num_pp; j++) { /**< Update positions of all birds for this process in parallel. */
            updateBirdPos(&proc_birds[j]);
        }

        MPI_Allgather(proc_birds, num_pp * BIRD_DOUBLES, MPI_DOUBLE, birds, num_pp * BIRD_DOUBLES, MPI_DOUBLE, MPI_COMM_WORLD); /**< Gather all birds' data. */

        #pragma omp parallel private(j)
        {
            #pragma omp for schedule(static)
            for (j = 0; j < num_pp; j++){ /**< Calculate angle effects for all birds for this process in parallel. */
                calculateAngleEffects(&proc_birds[j], birds, R);
            }

            #pragma omp barrier

            #pragma omp for schedule(static)
            for (j = 0; j < num_pp; j++) { /**< Update angles of all birds for this process in parallel. */
                updateBirdAngle(&proc_birds[j]);
                memcpy(&print_values[j * PRINT_DOUBLES], &proc_birds[j], PRINT_DOUBLES * sizeof(double)); /**< Copy the position and velocity doubles to a separate array for writing */
            }
        }

        MPI_Offset offset = (startnum + NUMBER * i) * PRINT_DOUBLES * sizeof(double); /**< Calculate where in file to write based on rank and timestep*/
        MPI_File_set_view(fh, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL); /**< Set Offset */
        MPI_File_write_all(fh, print_values, num_pp * PRINT_DOUBLES, MPI_DOUBLE, MPI_STATUS_IGNORE); /**< Write bytes to file */
        
    }
    if (rank == 0) {
        printf("Time Taken for %d Processes, %d Threads: %f", size, omp_get_max_threads(), omp_get_wtime() - startTime); /**< Print total simulation time on process 0. */
    }

    MPI_File_close(&fh);

    MPI_Finalize(); /**< Finalize MPI. */

    free(birds); /**< Free memory allocated for all birds. */
    free(proc_birds); /**< Free memory allocated for birds of this process. */

    return 0; /**< Return 0 to indicate successful completion. */
}