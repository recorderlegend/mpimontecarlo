#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

double f(double x);
double montecarlo(double a, double b, int iterations);
void input(int rank, int world, double *pointer_a, double *pointer_b, int *pointer_iterations);
void mpi_type(double *pointer_a, double *pointer_b, int *pointer_iterations, MPI_Datatype *mpi_input);
int main(int argc, char **argv)
{
    int rank, world;
    int n, local_n; // number of iterations, number of iterations per process
    double local_int, total_int;
    double a, b, dx, local_a, local_b; // start and end of integral             // change in x
    double startwtime = 0.0, endwtime;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    input(rank, world, &a, &b, &n);

    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }

    dx = (b - a) / n;
    local_n = n / world; // amount of iterations each process(world) handles

    local_a = a + rank * local_n * dx;                                                             // start of integral + (rank# * number of iterations * change in x)
    local_b = local_a + local_n * dx;                                                              // locala +(number of iterations * change in x)
    local_int = montecarlo(local_a, local_b, local_n);                                             // run monte carlo
    printf("n: %d | a: %f | b: %f | montecarlo val: %f \n", local_n, local_a, local_b, local_int); // debug print

    MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // end timer
        endwtime = MPI_Wtime();
        printf("wall clock: %lf\n", endwtime - startwtime);
    }
    MPI_Finalize();
}

double f(double x)
{
    return pow(x, 3);
}

double montecarlo(double a, double b, int iterations)
{
    time_t t;
    srand((unsigned)time(&t));
    double random, fVal;
    double sum = 0;

    int curIt = 0;

    while (curIt < iterations - 1)
    {

        // Select a random number within the limits of integration
        float ran = (float)(rand());
        random = a + (ran / RAND_MAX) * (b - a);

        // Sample the function's values
        fVal = f(random);

        // Add the f(x) value to the running sum
        sum += fVal;

        curIt++;
    }

    double estimate = (b - a) * sum / iterations;

    return estimate;
}

void mpi_type(double *pointer_a, double *pointer_b, int *pointer_iterations, MPI_Datatype *mpi_input)
{
    MPI_Aint address_a, address_b, addresss_iterations;
    MPI_Aint displacement_arr[3] = {0};
    int blocklengths_arr[3] = {1, 1, 1};
    MPI_Datatype types_arr[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Get_address(pointer_a, &address_a);
    MPI_Get_address(pointer_b, &address_b);
    MPI_Get_address(pointer_iterations, &addresss_iterations);
    displacement_arr[1] = address_b - address_a;
    displacement_arr[2] = addresss_iterations - address_a;
    MPI_Type_create_struct(3, blocklengths_arr, displacement_arr, types_arr, mpi_input);
    MPI_Type_commit(mpi_input);
}

void input(int rank, int world, double *pointer_a, double *pointer_b, int *pointer_iterations)
{
    MPI_Datatype mpi_input;
    mpi_type(pointer_a, pointer_b, pointer_iterations, &mpi_input);
    if (rank == 0)
    {
        printf("Enter the starting point,ending point, and number of iterations\n");
        scanf("%lf %lf %d", pointer_a, pointer_b, pointer_iterations);
    }
    MPI_Bcast(pointer_a, 1, mpi_input, 0, MPI_COMM_WORLD);
    MPI_Type_free(&mpi_input);
}