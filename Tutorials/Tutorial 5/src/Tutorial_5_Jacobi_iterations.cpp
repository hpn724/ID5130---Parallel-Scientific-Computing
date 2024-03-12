/*
Go to README.md and read the setting up for the tutorial

1. Develop a sequential program (in C/C++ or FORTRAN) for solving the system of linear equations described
above using Jacobi iteration method.

*/


#include<iostream>
#include<omp.h>
#include<math.h>
#include<stdio.h>

using namespace std;


// Function to initialize the matrix

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
double* vector_diff(double* x1, double* x2, int n)
{
    double* diff = new double[n];
    for(int i=0;i<n;i++)
    {
        diff[i]=x1[i]-x2[i];
    }
    return diff;
}

//Function to calculate the norm of a vector
template<typename T> double vector_norm(T x, int n)
{
    int sum=0;
    for(int i=0;i<n;i++)
    {
        sum+=pow(x[i],2);
    }
    return sqrt(sum);
}

//Jacobi iteration function
double* jacobi_iteration(double* x, double** a, double* b, int n, double threshold)
{
    double* xk=x; 
    double* xk1=initialize_sols_vector(n, 2);
    double* temp_vector;
    int counter=0;

    double current_error = vector_norm(vector_diff(xk1, xk, n),n)/vector_norm(xk, n); 
    
    while( abs(current_error) > threshold)
    {
        current_error = vector_norm(vector_diff(xk1, xk, n),n)/vector_norm(xk, n); 
        cout<<endl<<"Performing Jacobi Iteration for iteration number : "<<++counter<< ". Current error value : "<< current_error<<endl;
        temp_vector=xk1;
        for(int i=0;i<n;i++)
        {
            double temp_sum=0.0;

            for(int j=0;j<n;j++)
            {
                if(j==i)
                {
                    continue;
                }
                temp_sum+=a[i][j]*xk[j];
            }
            xk1[i]=(1/a[i][i])*(b[i]-temp_sum);
        }
        cout<<"xk at iteration number "<<counter<<" : "<<endl;
        show_vector(xk,n);
        cout<<"xk1 at iteration number "<<counter<<" : "<<endl;
        show_vector(xk1,n);

        xk=temp_vector;
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
    double epsilon;
    cin>>epsilon;

    double* solution=jacobi_iteration(x, a, b, n, epsilon);

    cout<<endl<<"The final solution for x : ";
    show_vector(solution, n);

    return 0;
}