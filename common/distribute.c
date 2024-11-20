/******************************************************************************
 *
 * File:        distribute.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Distributes indices in the interval [i1, i2] between slots in
 *              the interval [0, nslot_used - 1]. The calling slot has ID
 *              `myrank'. The results are stored in 6 global variables, with
 *              the following relations between them:
 *                my_number_of_iterations = my_last_iteration 
 *                                                      - my_first_iteration + 1
 *                my_number_of_iterations = number_of_iterations[rank]
 *                my_first_iteration = first_iteration[rank]
 *                my_last_iteration = last_iteration[rank]
 *              This code is based on distribute.[ch] in EnKF-C package.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>
#if defined(MPI)
#include <mpi.h>
#endif
#include "distribute.h"

int nslot_allocated = 0;
int my_number_of_iterations = -1;
int my_first_iteration = -1;
int my_last_iteration = -1;
int* number_of_iterations = NULL;
int* first_iteration = NULL;
int* last_iteration = NULL;

/** Distributes indices in the interval [i1, i2] between slots in the interval
 ** [0, nslot - 1]. The results are written to 6 global variables.
 * @param i1 Start of the interval
 * @param i2 End of the interval
 * @param nslot_used Number of slots (CPUs) to be be used
 * @param nslot_total Total number of slots (CPUs)
 * @param myrank ID of the process
 */
void distribute_iterations(int i1, int i2, int nslot_used, int nslot_total, int myrank)
{
    int niter = i2 - i1 + 1;
    int npp, i, j;

    assert(i2 >= i1);
    assert(nslot_used > 0 && nslot_used <= nslot_total);

#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (number_of_iterations == NULL) {
        number_of_iterations = malloc(nslot_total * sizeof(int));
        first_iteration = malloc(nslot_total * sizeof(int));
        last_iteration = malloc(nslot_total * sizeof(int));
        nslot_allocated = nslot_total;
    } else if (nslot_used > nslot_allocated) {
        number_of_iterations = realloc(number_of_iterations, nslot_total * sizeof(int));
        first_iteration = realloc(first_iteration, nslot_total * sizeof(int));
        last_iteration = realloc(last_iteration, nslot_total * sizeof(int));
        nslot_allocated = nslot_total;
    }

    npp = niter / nslot_used;
    j = niter % nslot_used;
    for (i = 0; i < j; ++i)
        number_of_iterations[i] = npp + 1;
    for (i = j; i < nslot_used; ++i)
        number_of_iterations[i] = npp;
    for (i = nslot_used; i < nslot_total; ++i)
        number_of_iterations[i] = 0;

    first_iteration[0] = i1;
    last_iteration[0] = i1 + number_of_iterations[0] - 1;
    for (i = 1; i < nslot_total; ++i) {
        first_iteration[i] = last_iteration[i - 1] + 1;
        last_iteration[i] = first_iteration[i] + number_of_iterations[i] - 1;
    }

    my_first_iteration = first_iteration[myrank];
    my_last_iteration = last_iteration[myrank];
    my_number_of_iterations = number_of_iterations[myrank];

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/**
 */
void distribute_free(void)
{
    if (number_of_iterations == NULL)
        return;
    free(number_of_iterations);
    number_of_iterations = NULL;
    free(first_iteration);
    first_iteration = NULL;
    free(last_iteration);
    last_iteration = NULL;
    my_number_of_iterations = -1;
    my_first_iteration = -1;
    my_last_iteration = -1;
}
