#include <cmath>
#include<iostream>
#include<math.h>
#include<omp.h>
#include<fstream>

using namespace std;


//two variable function defined
double func(double x, double y)
{
    return 2*(2-pow(x,2)-pow(y,2));
}

//double derivative of the defined two variable function
double d2f(double x, double y)
{
    return (pow(x,2)-1)*(pow(y,2)-1);
}

//function to show a matrix
void show_matrix(double** A, int m, int n) 
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}



//function to return the difference of two matrices
double** matrix_diff(double** A, double ** B, int m, int n)
{
    double **C=new double*[m] ; 
    
    for(int i=0;i<m;i++)
    {
        
        C[i]=new double[n] {};
        
        for(int j=0;j<n;j++)
        {
            C[i][j]=A[i][j]-B[i][j];
        }
    }
    return C;
}

//function to get the number of points in the grid 
int get_num_point_x(double start, double end, double grid_size)
{
    return int((end-start)/grid_size)+1;
}

//function to initialize the square grid
void initialize_grid(double start, int n, double grid_size, double* &x,double* &y)
{
    x = new double[n];
    y = new double[n];
    
    for(int i=0;i<n;i++)
    {   
        x[i] = start + i*grid_size;
        y[i] = start + i*grid_size;    
        
    }    
}

//fucntion to initialize the solution matrices
void initialize_solution_matrices(int n, double* x,double* y,double** phi, double** actual_values, double** q, double (*func)(double,double), double (*d2f)(double,double), double guess_value)
{
    phi = new double*[n];
    actual_values = new double*[n];
    q = new double*[n];

    for(int i=0;i<n;i++)
    {
        
        phi[i] = new double[n];
        actual_values[i] = new double[n];
        q[i] = new double[n];
        for(int j=0;j<n;j++)
        {
            q[i][j] = func(x[i], y[j]);
            actual_values[i][j] = d2f(x[i],y[j]);
            phi[i][j] = guess_value;
            
        }
    }


}

//function to find L2-norm of a matrix
double norm(double** A, int n)
{
    double res=0.0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            res+= A[i][j] * A[i][j];
            
        }
    }
    return sqrt(res);
}

//function to find error in Gauss Seidel iterations
double error_function(int n,double** actual, double** solved)
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
    return sqrt(error) / norm(actual, n);
}

void get_plot_data(double*x, double* y1, double* y2,int n)
{
    std::ofstream dataFile("Poisson_GaussSeidel_serial_data.csv");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        if(i==0)
        {
            dataFile << "x,y values,y solved" << endl;
        }
        dataFile << x[i] << "," << y1[i] << "," << y2[i] << endl;
    }
    dataFile.close();
}

//function to perform the Gauss Seidel iteration and return the final value matrix
double** gauss_seidel_iteration(int n,double grid_size,double** phi,double** actual_values,double** q,long long int iter_limit,double err_thresh, double (*error_function)(int, double**,double**))
{
    int count = 0;
    double error=100.0;

    double** phi_1 = new double*[n];

    phi_1 = new double*[n];
    for (int i = 0; i < n; i++) 
    {
        phi_1[i] = new double[n] {};
    }
    

    while(fabs(error)>err_thresh && iter_limit>count)
    {
        for(int i=1;i<n-1;i++)
        {
            
            for(int j=1;j<n-1;j++)
            {
                phi_1[i][j] = 0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1] + pow(grid_size,2)*q[i][j]);
            }
        }
        
        
        error = error_function(n,actual_values,phi_1);
        

        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                phi[i][j]=phi_1[i][j];
            }
        }

        count++;


        cout<<"Error value : "<< error << " at "<<"iteration no : "<<count<<endl;
        
    }

    if(abs(error)<=err_thresh)
    {
        cout<<"Convergence achieved at "<<count<<" number of iterations !!"<<endl;
    }
    
    else if(abs(error)>err_thresh && count > iter_limit)
    {
        cout<<"Iteration limit reached."<<endl;
    }

    return phi_1;
}


//function to obtain plot data



int main()
{

    double grid_size = 0.1;
    double start = -1;
    double end = 1;
    int n = static_cast<int>(((end-start)/grid_size))+1;

    double *x,*y;
    
    initialize_grid(start, n, grid_size,  x,  y);

    double **phi,  **actual_values, **q = new double*[n];
    
    double guess_value=0.0;


    initialize_solution_matrices(n,x,y,phi,actual_values,q,&func,&d2f,guess_value);

    double err_thresh = 0.01;

    long long int iter_limit=10000000;

    
    double** phi_1 = gauss_seidel_iteration(n,grid_size,phi,actual_values,q,iter_limit,err_thresh,&error_function);

    
    cout<<"Final value :"<<endl;
    show_matrix(phi_1, n,n);
    cout<<"Actual value :"<<endl;
    show_matrix(actual_values, n,n);


    // y =0.5 at i = 15

    double* ploty=new double[n];
    double* ploty_actual=new double[n];

    for(int i=0;i<n;i++)
    {
        ploty[i]=phi_1[i][15];
        ploty_actual[i]=actual_values[i][15];
    }

    get_plot_data(x,ploty_actual,ploty,n);

    for(int i=0;i<n;i++)
    {
        delete[] phi[i];
        delete[] phi_1[i];
        delete[] actual_values[i];
        delete[] q[i];
    }

    delete[] x;
    delete[] y;
    delete[] phi;
    delete[] phi_1;
    delete[] actual_values;
    delete[] q;

    return 0;
}