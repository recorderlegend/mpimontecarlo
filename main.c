#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double f(double x);
double montecarlo(double a, double b, int iterations);

int main(int argc, char **argv)
{
    int rank, world;
    int n, local_n;          // number of iterations, number of iterations per process
    n = 200;
    double local_int, total_int;
    double a = 1.0, b = 5.0; // start and end of integral
    double dx;               // change in x
    double local_a, local_b;
    int source;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    dx = (b-a)/n;
    local_n = n/world; //amount of iterations each process(world) handles
    local_a = a + rank*local_n*dx; // start of integral + (rank# * number of iterations * change in x)
    local_b = local_a + local_n*dx; //locala +(number of iterations * change in x)
    local_int = montecarlo(local_a, local_b, local_n); // run monte carlo
    printf("n: %d | a: %f | b: %f | montecarlo val: %f \n", local_n, local_a, local_b, local_int); // debug print
    if (rank != 0)
    {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        total_int = local_int;
        for (source = 1; source < world; source++)
        {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
        }
        printf("total integral: %f\n", total_int);
    }
    MPI_Finalize();

}

double f(double x)
{
	return pow(x,4) * exp(-x);
}

double montecarlo(double a, double b, int iterations){
    time_t t;
    srand((unsigned)time(&t));
    double totalSum = 0;
    double randNum, functionVal;

    int iter = 0;

    while (iter < iterations - 1)
    {

        // Select a random number within the limits of integration
        float ran = (float)(rand());
        randNum = a + (ran / RAND_MAX) * (b - a);

        // Sample the function's values
        functionVal = f(randNum);

        // Add the f(x) value to the running sum
        totalSum += functionVal;

        iter++;
    }

    double estimate = (b - a) * totalSum / iterations;

    return estimate;
}