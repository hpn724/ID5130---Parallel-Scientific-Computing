/* MPI parallel version of trapezoidal rule */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

#define pi 3.14159265358

double func(double x)
{
    return sin(x)/(2*pow(x,3));
}

double trapezoidal_rule(double la, double lb, double ln, double h)
{
    double total;
    double x;
    int i;

    total = (func(la) + func(lb))/2.0;
    for(i = 1; i <= ln-1; i++) /* sharing the work, use only local_n */
        {
            x = la + i*h;
            total += func(x);
        }
    total = total * h;

    return total;			/* total for each thread, private */
}


int main(int argc, char* argv[])
{
    double a, b, final_result, la, lb, lsum, h;
    int myid, nprocs, proc;
    int n, ln;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);	    /* myrank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);     /* size of the communicator */
    MPI_Status status;

    n = 256;			/* number of trapezoids.. */
    a = 1;
    b = pi;			    /* hard-coded.. */
    final_result = 0.0;

    h = (b-a)/n;
    ln = n/nprocs; 		/* nprocs evenly divides number of trapezoids */

    la = a + myid*ln*h;
    lb = la + ln*h;
    lsum = trapezoidal_rule(la, lb, ln, h); /* every process calls this function... */

    if (myid != 0)
    {
        MPI_Send(&lsum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    else				/* process 0 */
    {
        final_result = lsum;
        for (proc = 1; proc < nprocs; proc++)
        {
            MPI_Recv(&lsum, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);
            final_result += lsum;
        }

        printf("\n The area calculated under the curve sin(x)/(2*x^3) between 0 to PI using %d processors is equal to %lf \n\n", nprocs ,final_result);
    }
    
    MPI_Finalize();

    return 0; 
}
