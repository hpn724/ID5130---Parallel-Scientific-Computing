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
