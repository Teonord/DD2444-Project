#include "proj_common.h"
#include <stdio.h>
#include <omp.h>
#include <string.h>

/**
 * @brief Struct to represent the link between a bird and its box
 *
 * This struct connects a bird struct to its assigned box number
 */
struct Box
{
    struct Bird *bird;
    int boxnum;
};

/**
 * @brief function to compare two birds based on box number 
 * 
 * @param a Birdbox 1 to compare
 * @param b Birdbox 2 to compare
 * @return Difference in box number between the two bird boxes. 
 */
int compare(const void *a, const void *b) {
    struct Box *pa = (struct Box *)a;
    struct Box *pb = (struct Box *)b;
    return pa->boxnum - pb->boxnum;
}

/**
 * @brief Calculates the effects of neighboring birds on the current bird's angle.
 *
 * This function calculates the effects of neighboring birds within a certain sphere radius (R)
 * on the angle of the current bird. It updates the sum components presented in
 * the paper "Consensus of the 3-Dimensional Vicsek Model" by Liu Zhixin
 *
 * @param b Pointer to the bird struct to update its angle effects.
 * @param birds Array of all birds in the simulation.
 * @param R Pre-squared radius within which neighboring birds are considered.
 */
void boxCalculateAngleEffects(struct Bird *b, struct Bird **birds, double R, int start, int end) {
    struct Bird *nb; /**< Pointer to a neighboring bird. */
    for (int k = start; k < end; k++) {
        nb = birds[k];
        if (pow(nb->x - b->x, 2) + pow(nb->y - b->y, 2) + pow(nb->z - b->z, 2) < R) /**< Check if bird is within squared spherical radius R. */
        {
            b->sint += sin(nb->theta); /**< Update sum of sin theta. */
            b->costcose += cos(nb->theta) * cos(nb->eta); /**< Update sum of cos theta and eta. */
            b->costsine += cos(nb->theta) * sin(nb->eta); /**< Update sum of cos theta and sin eta. */
        }
    }
}

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
    
    int row = L / R_INIT; /**< Amount of boxes are in one row */
    double box_size = (double)L / (double)row; /** Side size of each box */
    int side = row * row; /**< Amount of boxes in 2D */
    int amount_boxes = side * row; /**< Total amount of boxes in 3D */
    int box_pos[amount_boxes]; /**< Amount of particles in each box */

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
    struct Bird **new_birds = calloc(NUMBER, __SIZEOF_POINTER__); /**< Allocate memory for all birds when sorted. */
    struct Bird *proc_birds = calloc(num_pp, sizeof(struct Bird)); /**< Allocate memory for birds of this process. */
    double *print_values = malloc(num_pp * PRINT_DOUBLES * sizeof(double)); /**< Allocate memory for values to save from this process. */

    struct Box *boxpairs = calloc(NUMBER, sizeof(struct Box));

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

        #pragma omp parallel for schedule(static) private(j)
        for (j = 0; j < amount_boxes; j++) box_pos[j] = 0;

        for (j = 0; j < NUMBER; j++){ /**< Calculate box number for all birds*/
            struct Bird *b = &birds[j];
            int boxnum = floor(b->x / box_size) + floor(b->y / box_size) * row + floor(b->z / box_size) * side;
            boxpairs[j].boxnum = boxnum;
            boxpairs[j].bird = b;
            if (boxnum + 1 != amount_boxes) { /**< Increment amounts of particles in each box. If pos is (L, L, L) it does not matter. */
                box_pos[boxnum + 1] += 1;
            }
        }

        for (int j = 0; j < amount_boxes - 1; j++) { /**< Makes sure each index shows the amount of particles before this box starts. */
            box_pos[j + 1] += box_pos[j];
        }

        qsort(boxpairs, NUMBER, sizeof(struct Box), compare); /**< Sort all birds to be in the order of box numbers */

        #pragma omp parallel private(j)
        {
            #pragma omp for schedule(static)
            for (j = 0; j < NUMBER; j++){ /**< Get list of birds in order of which box they are in */
                new_birds[j] = boxpairs[j].bird;
            }

            #pragma omp barrier
            
            int x, y, z;
            #pragma omp for schedule(static) private(x, y, z)
            for (j = 0; j < num_pp; j++){ /**< Calculate angle effects for all birds for this process in parallel. */
                struct Bird *b = &proc_birds[j];
                int cx = floor(b->x / box_size);
                if (cx == row) cx--;
                int cy = floor(b->y / box_size);
                if (cy == row) cy--;
                int cz = floor(b->z / box_size);
                if (cz == row) cz--;
                for (x = -1; x <= 1; x++) {
                    for (y = -1; y <= 1; y++) {
                        for (z = -1; z <= 1; z++) {
                            int nb_box = modd(cx + x, row) + modd(cy + y, row) * row + modd(cz + z, row) * side; /**< Calculate box index */
                            if (nb_box + 1 != amount_boxes) boxCalculateAngleEffects(b, new_birds, R, box_pos[nb_box], box_pos[nb_box + 1]); /**< Calculates distance to all birds in box, unless it is the last box */
                            else boxCalculateAngleEffects(b, new_birds, R, box_pos[nb_box], NUMBER); /**< If last box, just check rest of particles. */
                        }
                    }
                }
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
        MPI_File_write_all(fh, print_values, num_pp * PRINT_DOUBLES, MPI_DOUBLE, MPI_STATUS_IGNORE);  /**< Write bytes to file */
        
    }
    if (rank == 0) {
        printf("Time Taken for %d Processes, %d Threads: %f", size, omp_get_max_threads(), omp_get_wtime() - startTime); /**< Print total simulation time on process 0. */
    }

    MPI_File_close(&fh);

    MPI_Finalize(); /**< Finalize MPI. */

    free(birds); /**< Free memory allocated for all birds. */
    free(proc_birds); /**< Free memory allocated for birds of this process. */
    free(new_birds);
    free(boxpairs);

    return 0; /**< Return 0 to indicate successful completion. */
}