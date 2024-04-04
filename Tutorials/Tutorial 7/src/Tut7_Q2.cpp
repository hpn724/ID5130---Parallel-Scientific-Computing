#include<iostream>
#include<string.h>
#include<stdio.h>
#include "mpi.h"

using namespace std;

int main(int argc, char** argv)
{
    int i, myid, size, tag=100;
    int message_send, message_recv;
    MPI_Status status;
    int buffer = 1;


    //Initialize MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    if(myid!=0)
    {

        message_send = myid*(8-myid);                                   //number sending from different process
        MPI_Send(&message_send,buffer,MPI_INT,0,tag,MPI_COMM_WORLD);    //send message
    }

    else
    {
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"\t \t Running Tutorial 7 Question 2"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        int sum = 0;
        for(i=1; i<size; i++)
        {
            MPI_Recv(&message_recv,buffer,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
            printf("Number recieved from process id : %d = %d \n",i, message_recv);
            sum += message_recv;
        }
        printf("The final sum = %d \n",sum);
        
    }

    MPI_Finalize();

    return 0;
}