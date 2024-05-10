#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include "operations_omp.h"



using namespace std;


vector<vector<double>> glorot_init(int nin, int nout) 
{
    double sd = sqrt(6.0 / (nin + nout));
    vector<vector<double>> result(nin, vector<double>(nout));
    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<double> dist(-sd, sd);
    for (int i = 0; i < nin; ++i) 
    {
        for (int j = 0; j < nout; ++j) 
        {
            result[i][j] = dist(gen);
        }
    }

    return result;
}


//Function to display vector elements
template<typename T>
void show_vector(vector<T> A)
{
    cout<<"[";
    for(int i=0;i<A.size();i++)
    {
        if(i!=A.size()-1)
        {
            cout<<A[i]<<", ";
        }
        else
        {
            cout<<A[i];
        }
    }
    cout<<"]"<<endl;

} 

//Function to display matrix elements
void show_matrix(vector<vector<double>> A)
{
    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[i].size();j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}

//Function to flatten a matrix to a 1D element
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& nested_vector) 
{
    std::vector<T> flattened_vector;

    //#pragma omp parallel for collapse(2) shared(nested_vector, flattened_vector) schedule(static)
    for (size_t i = 0; i < nested_vector.size(); ++i) 
    {
        for (size_t j = 0; j < nested_vector[0].size(); ++j) 
        {
            flattened_vector.push_back(nested_vector[i][j]);
        }
    }

    return flattened_vector;
}

//Function to reshape 1D element to a matrix

vector<vector<double>> reshape(const vector<double>& input_vector, int rows, int cols) 
{
    vector<vector<double>> reshaped_matrix(rows, vector<double>(cols));

    int i,j,index = 0;
    
    //#pragma omp parallel for collapse(2) shared(input_vector, reshaped_matrix) schedule(static) 
    for (i = 0; i < rows; ++i) 
    {
        for (j = 0; j < cols; ++j) 
        {
            reshaped_matrix[i][j] = input_vector[index++];
            
        }
    }

    return reshaped_matrix;
}

//Function to find norm of a vector
double norm(vector<double> A)
{
    double sum=0;
    int i;
    #pragma omp parallel for reduction(+ : sum)
    for(i=0;i<A.size();i++)
    {
        sum+=pow(A[i],2);
    }
    return sqrt(sum);
}

//Function to find the difference between two vectors
vector<double> vector_diff(vector<double> A, vector<double> B)
{
    vector<double> diff(A.size());

    #pragma omp parallel for
    for (int i = 0; i < A.size(); i++) 
    {
        diff[i] = A[i] - B[i];
    }

    return diff;
}

//Function to find the normalized difference between two vectors
double norm_diff(vector<double> A, vector<double> B)
{
    return norm(vector_diff(A,B))/(norm(A)+norm(B));
}


//Function to calculate transpose of a matrix
vector<vector<double>> transpose(const vector<vector<double>>& matrix)
{
    vector<vector<double>> result(matrix[0].size(), vector<double>(matrix.size(), 0.0));

    int i,j;
    #pragma omp parallel for collapse(2)
    for (i = 0; i < matrix.size(); ++i) 
    {
        for (j = 0; j < matrix[0].size(); ++j) 
        {
            result[j][i] = matrix[i][j];
        }
    }
    
    
    
    return result;
    
}

//Function for the dot product of two vectors

double dotProduct(vector<double> A, vector<double> B)
{
    double sum=0;
    int i;
    #pragma omp parallel for shared(A,B) 
    for(i=0;i<A.size();i++)
    {
        #pragma omp atomic
        sum+=A[i]*B[i];
    }
    return sum;
}

//Function to calculate the mean of elements in a vector
template<typename T>
double mean(vector<T> A)
{
    double sum=0;
    int i;
    #pragma omp parallel for shared(A) reduction(+ : sum)
    for(int i=0;i<A.size();i++)
    {
        sum+=A[i];
    }
    return sum/A.size();
}

//Function to return a identity matrix
vector<vector<double>> identity_matrix(int n)
{
    vector<vector<double>> I(n,vector<double>(n,0));
    int i;
    #pragma omp parallel for
    for(i=0;i<n;i++)
    {
        I[i][i] = 1;
    }
    return I;
}

//Function to calculate the square root inverse of a diagonal matrix
vector<vector<double>> square_root_inverse_diag_matrix(vector<vector<double>> diag)
{    
    size_t n = diag.size();
    vector<vector<double>> inv(n, vector<double>(n, 0.0));
    int i;
    #pragma omp parallel for 
    for (i = 0; i < n; i++) 
    {
        
        inv[i][i] = 1.0 / sqrt(diag[i][i]);
        
    }
    return inv;
}


//Function to add two matrices
vector<vector<double>> matrix_sum(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t rows_B = B.size();
    size_t cols_B = B[0].size();

    // Determine the dimensions of the result matrix
    size_t rows_result = max(rows_A, rows_B);
    size_t cols_result = max(cols_A, cols_B);

    // Initialize the result matrix with zeros
    vector<vector<double>> result(rows_result, vector<double>(cols_result, 0.0));
    int i,j;
    #pragma omp parallel for collapse(2)
    for ( i = 0; i < rows_result; ++i) 
    {
        for ( j = 0; j < cols_result; ++j) 
        {
            if (i < rows_A && j < cols_A) 
            {
                // Add element from A
                #pragma omp atomic
                result[i][j] += A[i][j];
            }
            if (i < rows_B && j < cols_B) 
            {
                // Add element from B
                #pragma omp atomic
                result[i][j] += B[i][j];
            }
        }
    }
    

    return result;
}



//Function to multiply two matrices
vector<vector<double>> matrix_product(const vector<vector<double>> A, const vector<vector<double>> B)
{
    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t rows_B = B.size();
    size_t cols_B = B[0].size();

    if (cols_A != rows_B) {
        cerr << "Error: Incompatible dimensions for matrix multiplication." << endl;
        return {};
    }

    double temp ;
    vector<vector<double>> product(rows_A, vector<double>(cols_B, 0.0));
    int i,j,k;
    #pragma omp parallel for private(i,j,k) shared(rows_A,cols_B,cols_A,A,B,product) reduction(+ : temp)
    for (i = 0; i < rows_A; ++i) 
    {
        for (j = 0; j < cols_B; ++j) 
        {
            temp = 0;
            for (k = 0; k < cols_A; ++k) 
            {
                temp += A[i][k] * B[k][j];
            }
            product[i][j] = temp;
        }
    }


    return product;
}
