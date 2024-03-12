#include<iostream>
#include<vector>
#include<random>
#include<omp.h>

using namespace std;

vector<vector<int>> get_random_matrix(int n, int limit)
{
    vector<vector<int>> matrix ;
    for (int i = 0; i < n; ++i) 
    {
        vector<int> row;
        for (int j = 0; j < n; ++j) {
            row.push_back(rand()%limit);
        }
        matrix.push_back(row);
    }
    return matrix;
}

void show_matrix(vector<vector<int>> mat)
{
    for(int i=0;i<mat.size();i++)
    {
        for(int j=0;j<mat[0].size();j++)
        {
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    return ;
}

pair<int**, double> parallel_add(vector<vector<int>> A, vector<vector<int>> B, int threads)
{
    int** C=new int[A.size()];
    int i=0,j=0;
    double t;
    t=omp_get_wtime();
    #pragma omp for default(none) collapse(2) shared(A,B,C) private(i,j) num_threads(threads)
    for(i=0;i<A.size();i++)
    {
        C[i]=new int[A[0].size()];
        for(j=0;j<A[0].size();j++)
        {
            C[i][j]=A[i][j]+B[i][j];
        }
    }
    t=omp_get_wtime()-t;
    pair<int**, double> ans(C,t);
    return ans;
}


int main() 
{
    int n;
    cout<<"Enter size of matrix : ";
    cin>>n;
    int limit;
    cout<<"Enter element upper limit of the matrix : ";
    cin>>limit;
    vector<vector<int>> A=get_random_matrix(n,limit);
    vector<vector<int>> B=get_random_matrix(n, limit);
    pair<int**, double> C=parallel_add(A, B, 2);

    show_matrix(A);
    show_matrix(B);

    cout<<"The addition result , C :"<<endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout<<C.first[i][j]<<" ";
        }
        cout<<endl;
    }

    cout<<"The time taken, t : "<<C.second<<endl;


    return 0;
}