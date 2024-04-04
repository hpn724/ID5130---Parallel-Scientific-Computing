#include<iostream>
#include<string.h>
#include<stdio.h>
#include "mpi.h"

using namespace std;

int main(int argc, char** argv)
{
    int i, myid, size, tag=100;
    char message_send[50], message_recv[50];
    MPI_Status status;


    //Initialize MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);



    if(myid == 0)
    {

        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"\t \t Running Tutorial 7 Question 1"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        sprintf(message_send, "Calling from process id : %d \n",myid);
        for(int i=1; i<size; i++)
        {
            MPI_Send(message_send,50,MPI_CHAR,i,tag,MPI_COMM_WORLD);
        }
    }
        

    else
    {
        
        MPI_Recv(message_recv,50,MPI_CHAR,0,tag,MPI_COMM_WORLD,&status);
        printf("Message received to process %d - %s ", myid, message_recv);
    }

    MPI_Finalize();

    return 0;
}