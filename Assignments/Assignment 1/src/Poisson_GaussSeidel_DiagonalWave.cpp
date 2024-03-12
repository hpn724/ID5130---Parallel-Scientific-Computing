#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>

using namespace std;

// Function to calculate the two-variable function
double func(double x, double y) 
{
    return 2 * (2 - pow(x, 2) - pow(y, 2));
}

// Function to calculate the double derivative of the function
double d2f(double x, double y) 
{
    return (pow(x, 2) - 1) * (pow(y, 2) - 1);
}

void show_matrix(double** A, int m, int n) 
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

// Function to initialize the square grid
void initialize_grid(double start, int n, double grid_size, double*& x, double*& y) 
{
    x = new double[n];
    y = new double[n];
    #pragma omp parallel for
    for (int i = 0; i < n; i++) 
    {
        x[i] = start + i * grid_size;
        y[i] = start + i * grid_size;
    }
}

// Function to initialize the solution matrices
void initialize_solution_matrices(int n, double* x, double* y, double**& phi, double**& actual_values, double**& q, double (*func)(double, double), double (*d2f)(double, double), double guess_value, int thread) 
{
    phi = new double*[n];
    actual_values = new double*[n];
    q = new double*[n];
    #pragma omp parallel for shared(phi, actual_values, q) num_threads(thread)
    for (int i = 0; i < n; i++) 
    {
        phi[i] = new double[n];
        actual_values[i] = new double[n];
        q[i] = new double[n];
        for (int j = 0; j < n; j++) 
        {
            q[i][j] = func(x[i], y[j]);
            actual_values[i][j] = d2f(x[i], y[j]);
            phi[i][j] = guess_value;
        }
    }
}


// Function to calculate the L2-norm of a matrix
double norm(double** A, int n, int thread) 
{
    double res = 0.0;
    #pragma omp parallel for num_threads(thread) reduction(+:res)
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            res += A[i][j] * A[i][j];
        }
    }
    return sqrt(res);
}

// Function to calculate the error in Gauss-Seidel iterations
double error_function(int n, double** actual, double** solved, int thread) 
{
    double error = 0.0;
    #pragma omp parallel for num_threads(thread) reduction(+:error)
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            error += pow(actual[i][j] - solved[i][j], 2);
        }
    }
    return sqrt(error) / norm(actual, n,thread);
}

// Function to obtain plot data
void get_plot_data(double* x, double* y1, double* y2, int n, double grid_size, int thread) 
{
    string file_name = "Poisson_GaussSeidel_DiagonalWave_data_grid_" + to_string(grid_size) + "_threads_" + to_string(thread) + ".csv";
    ofstream dataFile(file_name);
    if (!dataFile.is_open()) 
    {
        cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    dataFile << "x,y values,y solved" << endl;
    for (int i = 0; i < n; i++) 
    {
        dataFile << x[i] << "," << y1[i] << "," << y2[i] << endl;
    }
    dataFile.close();
}

// Function to obtain grid time data
void get_grid_time_data(double* x, double* y) 
{
    string file_name = "Poisson_GaussSeidel_DiagonalWave_data_grid_time_data.csv";
    ofstream dataFile(file_name);
    if (!dataFile.is_open()) 
    {
        cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    dataFile << "grid size,time taken" << endl;
    for (int i = 0; i < 3; i++) 
    {
        dataFile << x[i] << "," << y[i] << endl;
    }
    dataFile.close();
}

// Function to obtain grid time data
void get_thread_time_data(int* x, double* y) 
{
    string file_name = "Poisson_GaussSeidel_DiagonalWave_data_thread_time_data.csv";
    ofstream dataFile(file_name);
    if (!dataFile.is_open()) 
    {
        cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    dataFile << "grid size,time taken" << endl;
    for (int i = 0; i < 3; i++) 
    {
        dataFile << x[i] << "," << y[i] << endl;
    }
    dataFile.close();
}


// Function to perform the Gauss-Seidel iteration
double gauss_seidel_diag_wave_iteration(int n, double*x, double grid_size, double**& phi, double** actual_values, double** q, long long int iter_limit, double err_thresh, int thread, double* timer, double (*error_function)(int, double**,double**,int)) 
{
    double start_time = omp_get_wtime();
    double error=100.0;
    int count = 0;
    double** phi_1 = new double*[n];
    for (int i = 0; i < n; i++) 
    {
        phi_1[i] = new double[n] 
        {};
    }

    // Gauss-Seidel iteration loop
    while (count < iter_limit && fabs(error) > err_thresh) 
    {
        #pragma omp parallel for shared(phi, phi_1)
        for (int l = 1; l < 2 * n - 1; l++) 
        {
            int ibeg = max(1, l - n + 2);
            int iend = min(n - 2, l - 1);
            for (int i = ibeg; i <= iend; i++) 
            {
                int j = l - i;
                phi_1[i][j] = 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] + grid_size * grid_size * q[i][j]);
            }
        }
        error = error_function(n,actual_values,phi_1,thread);
        #pragma omp parallel for
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                phi[i][j] = phi_1[i][j];
            }
        }
        count++;
    }



    double end_time = omp_get_wtime();
    *timer = end_time - start_time;

    if(abs(error)<=err_thresh)
    {
        printf("Convergence achieved for grid size %f with %d number of threads in %f seconds !!\n", grid_size, thread, *timer);

    }
    else if(abs(error)>err_thresh && count > iter_limit)
    {
        printf("Iteration limit reached \n");
    }
    else
    {
        printf("Convergence not achieved for grid size %f with %d number of threads in %f seconds :( \n", grid_size, thread, *timer);
    }


    if(grid_size==0.1 && thread ==2)
    {
        cout<<"Final value :"<<endl;
        show_matrix(phi_1, n,n);
        cout<<"Actual value :"<<endl;
        show_matrix(actual_values, n,n);

        double* ploty=new double[n];
        double* ploty_actual=new double[n];

        for(int i=0;i<n;i++)
        {
            ploty[i]=phi_1[i][15];
            ploty_actual[i]=actual_values[i][15];
        }
        get_plot_data(x,ploty_actual,ploty,n,grid_size,thread);
    
    }

    // Clean up memory
    for (int i = 0; i < n; i++) 
    {
        delete[] phi_1[i];
    }

    delete[] phi_1;

    return *timer;
}

double diagonal_wave_method_timer(double grid_size, int thread) 
{
    double start = -1;
    double end = 1;
    double* x;
    double* y;
    int n = static_cast<int>(((end-start)/grid_size))+1;
    initialize_grid(start, n, grid_size, x, y);

    double** phi;
    double** actual_values;
    double** q;
    initialize_solution_matrices(n, x, y, phi, actual_values, q, &func, &d2f, 0.0, thread);

    double timer;
    gauss_seidel_diag_wave_iteration(n, x,grid_size, phi, actual_values, q, 10000000, 0.01, thread, &timer,&error_function);

    // Cleanup
    delete[] x;
    delete[] y;
    for (int i = 0; i < n; i++) 
    {
        delete[] phi[i];
        delete[] actual_values[i];
        delete[] q[i];
    }
    delete[] phi;
    delete[] actual_values;
    delete[] q;

    return timer;
}

int main() 
{

    // Define grid sizes and thread counts
    double grid_sizes[] = {0.1, 0.01, 0.005};
    int num_threads[] = {2, 4, 8};
    double* grid_time = new double[3];


    //test for grid size=0.1 and thread=2
    double temp_time = diagonal_wave_method_timer(0.1, 2);


    // Iterate over grid sizes and thread counts
    for (int i=0;i<3 ; i++) 
    {
        
        double timer = diagonal_wave_method_timer(grid_sizes[i], num_threads[2]);
        grid_time[i]=timer;
        
    }

    get_grid_time_data(grid_sizes, grid_time);

    double* thread_time= new double[3];
    for (int i=0;i<3 ; i++) 
    {
        
        double timer = diagonal_wave_method_timer(grid_sizes[2], num_threads[i]);
        thread_time[i]=timer;
        
    }

    get_thread_time_data(num_threads, thread_time);


    return 0;
}