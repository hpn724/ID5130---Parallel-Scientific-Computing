#include<iostream>
#include<random>

using namespace std;

//Function to get the random values for matrix. Function takes the size and the maximum limit of the elements of matrix as arguments.
int** get_random_matrix(int m, int n, int limit)
{
    int** matrix;
    matrix=new int*[m];
    for(int i=0;i<m;i++)
    {
        matrix[i]=new int[n];
    }

    for (int i = 0; i < m; ++i) 
    {
        
        for (int j = 0; j < n; ++j) {
            matrix[i][j]=rand()%limit;
        }
    }
    return matrix;
}

//Function to get the transpose of the matrix
int** transpose(int** mat, int m, int n)
{
    int** trans;
    trans=new int*[n];
    for(int i=0;i<n;i++)
    {
        trans[i]=new int[n];
    }
    for(int i=0; i<m;i++)
    {
        
        for(int j=0;j<n;j++)
        {
            trans[i][j]=(mat[j][i]);
        }
    }
    return trans;
}

void show_matrix(int** mat, int m, int n)
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
    return ;
}

int** matrix_add(int** A, int m1, int n1, int** B, int m2, int n2)
{
    if(m1==m2 && n1==n2)
    {
        int** add;
        add = new int*[m1];
        for(int i=0;i<m1;i++)
        {
            add[i]=new int[n1];
        }

        for(int i=0;i<m1;i++)
        {
            
            for(int j=0;j<m2;j++)
            {
                add[i][j]=A[i][j]+B[i][j];
            }
        }
        return add;
    }
    else
    {
        cout<<"Matrix dimensions are not compatible for addition. Please check again "<<endl;
        return {};
    }
}


int main()
{
    int n;
    cout<<"Enter the square matrix size : ";
    cin>>n;
    
    int limit;
    cout<<"Enter the maximum allowed value for matrix element : ";
    cin>>limit;
    
    int** A=get_random_matrix(n, n, limit);
    int** B=get_random_matrix(n, n, limit);

    cout<<"A : "<<endl;
    show_matrix(A,n,n);

    int** sym=matrix_add(A, n, n, transpose(A,n,n), n, n);
    cout<<"A+transpose(A) : "<<endl;
    show_matrix(sym,n,n);
    
    return 0;
}