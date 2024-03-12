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

int** product(int** A, int m1,int n1, int** B,int m2, int n2)
{
    int** prod;
    prod=new int*[m1];
    for(int i=0;i<m1;i++)
    {
        prod[i]=new int[n2];
    }
    if(n1==m2)
    {
        for(int i=0;i<m1;i++)
        {
            for(int j=0;j<n2;j++)
            {
                int sum=0;
                for(int k=0;k<n1;k++)
                {
                    sum+=A[i][k]*B[k][j];
                }
                prod[i][j]=sum;
            }
        }
        return prod;
    }
    else
    {
        cout<<"The dimensions of the matrices specifies is not suitable for multiplication. Check the dimensions or try reversing the order"<<endl;
        return {};
    }
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


int main()
{
    int n;
    cout<<"Enter the square matrix size : ";
    cin>>n;
    
    int limit;
    cout<<"Enter the maximum allowed value for matrix element : ";
    cin>>limit;
    
    int** A=get_random_matrix(n,n,limit);
    int** B=get_random_matrix(n,n,limit);

    cout<<"A :"<<endl;
    show_matrix(A,n,n);


    cout<<"B :"<<endl;
    show_matrix(B,n,n);

    int** trans_prod_1=transpose(product(A,n,n,B,n,n),n,n);

    cout<<"transpose(A.B) : "<<endl;
    show_matrix(trans_prod_1,n,n);

    int** trans_prod_2=product(transpose(B,n,n),n,n, transpose(A,n,n),n,n);

    cout<<"transpose(B).transpose(A) : "<<endl;
    show_matrix(trans_prod_2,n,n);

    return 0;
}