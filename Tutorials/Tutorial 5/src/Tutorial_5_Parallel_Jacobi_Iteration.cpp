/*
Go to README.md and read the setting up for the tutorial

1. Develop a parallel program (in C/C++ or FORTRAN) for solving the system of linear equations described
above using Jacobi iteration method.

*/


#include<iostream>
#include<omp.h>
#include<math.h>
#include<stdio.h>
#include<algorithm>
#include<omp.h>


using namespace std;


// Function to initialize the matrix

double* matrix_vector_product(int n,double** a, double* v)
{
    double* result = new double[n]{};
    double temp;
    for(int i=0;i<n;i++)
    {
        temp = 0;
        for(int j=0;j<n;j++)
        {
            temp+= a[i][j]*v[j];
        }
        result[i]+=temp;
    }
    return result;
}

double** matrix_initializer(int m,int n)
{
    
    cout<<"Initialising a matrix."<<endl;
    double** mat;
    mat=new double*[m];
    for(int i=0;i<m;i++)
    {
        mat[i]=new double[n];
    }

    for(int i=1;i<=m;i++)
    {
        for(int j=1;j<=n;j++)
        {
            if (i==j)
            {
                mat[i-1][j-1]=i+j;
            }
            
            else if (i==1 && j==n)
            {
                mat[i-1][j-1]=1;
            }

            else if (i==n && j==1)
            {
                mat[i-1][j-1]=2*n-1;
            }

            else 
            {
                mat[i-1][j-1]=1.0/n;
            }
        }
    }
    return mat;
    
}

//Function to intialize the RHS side of the linear system of equation
double* RHS_initializer(int n)
{
    cout<<"Initialising RHS vector."<<endl;
    double* b = new double[n];
    for(int i=0;i<n;i++)
    {
        if(i==0)
        {
            b[i]=1;
        }
        else
        {
            b[i]=0;
        }
    }
    return b;
}

//Function to display a matrix 
void show_matrix(double** mat, int m, int n)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

//Function to initialize the vector x
double* initialize_sols_vector(int n,double value)
{
    cout<<"Initialising x vector."<<endl;
    double* x=new double[n];
    for(int i=0;i<n;i++)
    {
        x[i]=value;
    }
    return x;
}

//Function to display a vector
void show_vector(double* x, int n)
{
    cout<<"[";
    for(int i=0;i<n;i++)
    {
        if(i!=n-1)
        {
            cout<<x[i]<<", ";
        }
        else
        {
            cout<<x[i]<<"]";
        }
    }
    cout<<endl;
}

//Function to calculate the difference of two given vectors
double* vector_diff(double* x1, double* x2, int n, int threads)
{

    double* diff = new double[n];
    #pragma omp parallel for num_threads(threads)
    for(int i=0;i<n;i++)
    {
        diff[i]=x1[i]-x2[i];
    }
    return diff;
}

//Function to calculate the norm of a vector
double vector_norm(double* x, int n, int threads)
{
    double sum=0;
    int i;
    #pragma omp parallel for num_threads(threads) shared(i,x) reduction(+ : sum)
    for(i=0;i<n;i++)
    {
        sum+=pow(x[i],2);
    }
    return sqrt(sum);
}

//Jacobi iteration function
double* jacobi_iteration(double** a, double* b, int n, double init_value, double threshold, int threads)
{

    double time_taken = omp_get_wtime();
    
    double* xk=initialize_sols_vector(n,init_value); 
    double* xk1=new double[n];
    int counter=0;
    int iter_limit = pow(10,7);

    double current_error = 100; 
    
    int i,j;
    double temp_sum;

    #pragma omp parallel num_threads(threads) private(i,j,temp_sum) shared(xk1,xk,threads,a,b,current_error,counter,n,iter_limit,threshold) 
    {
        while( abs(current_error) > threshold && counter< iter_limit)
        {

            for(i=0;i<n;i++)
            {
                temp_sum=0.0;
                #pragma omp parallel for
                    for(j=0;j<n;j++)
                    {
                        if(j==i)
                        {
                            continue;
                        }
                        temp_sum+=a[i][j]*xk[j];
                    }

                xk1[i]=(1/a[i][i])*(b[i]-temp_sum);
            }

            current_error = vector_norm(vector_diff(xk1, xk, n, threads),n, threads)/vector_norm(xk, n, threads);

            cout<<endl<<"Iteration : "<<++counter<< "\t Current error value : "<< current_error<<endl;

            for(int i=0;i<n;i++)
            {
                xk[i]=xk1[i];
            }
        }
    }
    
    time_taken = omp_get_wtime() - time_taken;

    if(current_error<threshold)
    {
        cout<<"Solution has successfully converged within the threshold value!!"<<endl;
        cout<<"Time taken for convergence by "<<threads<<" number of threads : "<<time_taken<<endl;
    }

    else
    {
        cout<<"Solution failed to converge within threshold value for "<<counter<<" number of iterations. :("<<endl;
    }

    
    return xk1;
}

int main()
{
    cout<<endl;

    int n;
    cout<<"Enter the size of square matrix : ";
    cin>>n;
    cout<<endl;

    double** a=matrix_initializer(n,n);
    cout<<"The matrix a : "<<endl;
    show_matrix(a, n, n);

    double* b=RHS_initializer(n);
    cout<<"The vector b : "<<endl;
    show_vector(b, n);

    cout<<endl<<"Enter the value for initializing solution vector : ";
    double value;
    cin>>value;
    cout<<endl;
    
    double* x=initialize_sols_vector(n, value);
    cout<<"The initial solution vector x :"<<endl;
    show_vector(x, n);

    cout<<endl<<"Enter the threshold value : ";
    double epsilon = 0.0001;
    cin>>epsilon;


    // defining the number of threads to be used
    int threads[] = {2,4,8};

    double* solution = nullptr;
    for(auto thread : threads)
    {
        solution=jacobi_iteration(a, b, n,value, epsilon, thread);

        cout<<endl<<"The final solution for x : ";
        show_vector(solution, n);
        cout<<endl;


        cout<< "Verifying the result. The calculated value of b using solution obtained : "<<endl;
        show_vector(matrix_vector_product(n,a,solution),n);
        cout<<endl;
    }

    for(int i=0;i<n;i++)
    {
        delete[] a[i];
    }
    delete[] a;
    delete[] b;
    delete[] x;
    delete[] solution;

    return 0;
}