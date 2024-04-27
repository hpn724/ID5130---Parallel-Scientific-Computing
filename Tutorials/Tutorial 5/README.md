TUTORIAL 5

In this tutorial, you will write a serial and an OpenMP programs for the solution of Jacobi iteration method.
Consider the solution of system of linear equations, 𝐴𝑥 = 𝑏, where 𝐴 is the coefficient matrix, 𝑥 is the unknown
vector to be determined and 𝑏 is a known right-hand-side (RHS) vector. The matrix 𝐴 considered to be a dense and
diagonally dominant matrix given as follows:

𝐴 = [𝑎𝑖𝑗]𝑛×𝑛 and 1 ≤ 𝑖, 𝑗 ≤ 𝑛,

𝑎𝑖𝑗 =
            𝑖 + 𝑗, 𝑖 = 𝑗,
            1, 𝑖 = 1, 𝑗 = 𝑛,
            2𝑛 − 1, 𝑖 = 𝑛, 𝑗 = 1,
            1/𝑛, other cases,
 


The vector 𝑏 is given by 𝑏 = (1, 0, 0, ..., 0)𝑇 . The point Jacobi method uses the following equation to obtain an
improved solution, 𝑥𝑖 𝑘+1 from the guess/previous iteration value of 𝑥𝑖 𝑘 :



xk1[i] = (1/a[i][i])*(b[i] - sum[j=1][j=n](a[i][j]*xk[j]))

The iterations can begin from an initial guess of xk[i] and continue till the solution has converged to some tolerance (e) as per the following criterion

                        ||xk1 - xk||
                        ____________  <= e
                            ||xk||


1. Develop a sequential program (in C/C++ or FORTRAN) for solving the system of linear equations described above using Jacobi iteration method.

2. Develop an OpenMP version of the serial program that you developed earlier by suitably modifying the code. Convince yourself that the serial and the parallel programs that you developed are indeed giving correct result. Using the size of the system as n=10, 100 and 1000 and the number of threads as p = 2, 4 and 8 time your program. Do you notice any performance improvement using the parallel version of the program?