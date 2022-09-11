#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double f(double x);
double montecarlo(double a, double b, int iterations);

// Inspired by Pacheco's Parallel Processing algorithm for trapezoid integration
int main(int argc, char **argv)
{
    int rank, world;
    int n, local_n;          // number of iterations, number of iterations per process
    n = 10000000;
    double local_int, total_int;
    double a = 0.0, b = 5.0; // start and end of integral
    double dx;               // change in x
    double local_a, local_b;
    int source;
    double startwtime = 0.0, endwtime;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    if(rank == 0){
        startwtime = MPI_Wtime();
    }


    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    dx = (b-a)/n;
    local_n = n/world; //amount of iterations each process(world) handles
    
    
    local_a = a + rank*local_n*dx; // start of integral + (rank# * number of iterations * change in x)
    local_b = local_a + local_n*dx; //locala +(number of iterations * change in x)
    local_int = montecarlo(local_a, local_b, local_n); // run monte carlo
    printf("n: %d | a: %f | b: %f | montecarlo val: %f \n", local_n, local_a, local_b, local_int); // debug print

    
    // if (rank != 0)
    // {
    //     MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    // }
    // else
    // {
    //     total_int = local_int;
    //     for (source = 1; source < world; source++)
    //     {
    //         MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         total_int += local_int;
    //     }
    //     printf("total integral: %f\n", total_int);
    // }

    MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



    if(rank == 0){
        // end timer
        endwtime = MPI_Wtime();
        printf("wall clock: %lf\n", endwtime-startwtime);
    }
    MPI_Finalize();

}

double f(double x)
{
	return pow(x,3);
}

double montecarlo(double a, double b, int iterations)
{
    time_t t;
    srand((unsigned)time(&t));
    double sum = 0;
    double random, fval;
    int curit = 0;
    while (curit < iterations - 1)
    {
        // Select a random number within the limits of integration
        float ran = (float)(rand());
        random = a + (ran / RAND_MAX) * (b - a);
        // Sample the function's values
        fval = f(random);
        // Add the f(x) value to the running sum
        sum += fval;
        curit++;
    }
    double estimate = (b - a) * sum / iterations;
    return estimate;
}
