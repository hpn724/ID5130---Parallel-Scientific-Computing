#include<iostream>
#include<string.h>
#include<stdio.h>
#include "mpi.h"

using namespace std;

void show_vector(int n,int* arr)
{
    for(int i=0;i<n;i++)
    {
        if(i==0)
        {
            cout<<"["<<arr[i]<<", ";
        }
        else if(i==n-1)
        {
            cout<<arr[i]<<"]"<<endl;
        }
        else
        {
            cout<<arr[i]<<", ";
        }
    }
}

int main(int argc, char** argv)
{
    int i, j, myid, num_pros, tag=100;
    MPI_Status status;
    

    //Initialize MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_pros);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);





    if(myid == 0)
    {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"\t \t Running Tutorial 7 Question 3"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;

        
        int n = (num_pros-1)*4;
        int ln = n/(num_pros-1);

        int *local_arr = new int[ln]{};
        
    
        int *arr = new int[n]{};
        
        
        int *recv_back_arr = new int[ln]{};


        for(i=0;i<n;i++)
        {
            arr[i] = i;
        }
        
        printf("The initial array in process %d : ", myid);

        show_vector(n,arr);
        
        
        for(i=1;i<num_pros;i++)
        {
            for(j=0;j<ln;j++)
            {
                local_arr[j] = arr[(i-1)*ln+j];
            }
            printf("The message array send to process %d from process %d : ",i,myid);
            show_vector(ln, local_arr);
            MPI_Send(local_arr,ln,MPI_INT,i,tag,MPI_COMM_WORLD);

        }

        for(i=1;i<num_pros;i++)
        {
            MPI_Recv(recv_back_arr, n, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            
            for(j=0;j<n;j++)
            {
                arr[(i-1)*ln+j] = recv_back_arr[j];
            }
        }

        printf("The altered message array received back to process  %d : ",myid);
        show_vector(n,arr);


        delete[] arr;
        delete[] recv_back_arr;
        delete[] local_arr;

    }
        

    else
    {   
        int n = (num_pros-1)*4;
        int ln = n/(num_pros-1);
        int* recv_arr = new int[ln]{};

        MPI_Recv(recv_arr,ln,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
        
        for(int i=0;i<ln;i++)
        {
            recv_arr[i] += (myid-1)*recv_arr[i];
        }
        printf("The message array send to process %d from process %d : ", 0,myid);
        show_vector(ln,recv_arr);

        MPI_Send(recv_arr,ln,MPI_INT,0,tag,MPI_COMM_WORLD);
        delete[] recv_arr;
    }

    
    
    
    
    

    
    MPI_Finalize();

    return 0;
}