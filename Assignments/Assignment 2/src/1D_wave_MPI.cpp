#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<string.h>
#include "mpi.h"

using namespace std;

#define pi 3.14159265358

// function to show the vector (to verify the output if the o/p is correct or not)
void show_vector(int n, double* x)
{
    if(n==1)
    {
        cout<<"["<<x[0]<<"]"<<endl;
    }
    else
    {
        for(int i=0;i<n;i++)
        {
            if(i==0)
            {
                cout<<"["<<x[i];
            }
            else if(i==n-1)
            {
                cout<<", "<<x[i]<<"]"<<endl;
            }
            else
            {
                cout<<", "<<x[i];
            }
        }
    }
    
}

// Initial wave function
double initial_wave_function(double x,double left_pt, double mid_pt, double right_pt)
{
    if(x<=mid_pt && x>=left_pt)                             //condition for the sin(4*pi*x) value to be updated
    {
        return sin(4*pi*x);
    }
    else 
    {
        return 0.0;
    }


}

// Function to pass the inital wave points (t=0) as a list
double* initialize_wave(int n, double left_pt, double mid_pt, double right_pt, double dx)
{
    double* u0 = new double[n];
    double curr_x = left_pt;
    for(int i=0;i<n;i++)
    {        
        u0[i] = initial_wave_function(curr_x,left_pt, mid_pt,right_pt);
        curr_x+=dx;
    }
    return u0;
}

// Function to calculate the wave position at time t analytically
double* analytical_solution(int n,double local_start,double left_pt, double mid_pt, double right_pt, double dx, double c, double t)
{
    double* u = new double[n];
    double curr_x = local_start;
    for(int i=0;i<n;i++)
    {
        u[i] = initial_wave_function((curr_x-(c*t)),left_pt,mid_pt,right_pt);
        curr_x+=dx;
    }
    return u;
}

// Upwind scheme
double* upwind_scheme(int n,double* &u, double dx, double dt, double c, int proc_id, double* bwd_shdw)
{
    
    double* upd_u = new double[n];

    if(proc_id==0)
    {
        for(int i=1;i<n;i++)
        {
            upd_u[i] = u[i]-(c*dt/dx)*(u[i]-u[i-1]);
        }
    }

    else
    {
        double prev_proc_last_pt = bwd_shdw[(sizeof(bwd_shdw))/sizeof(bwd_shdw[0])];

        upd_u[0] = u[0]-(c*dt/dx)*(u[0]- prev_proc_last_pt);      //for points starting in between the value needs to be updated with the shadow points from preceding processor

        for(int i=1;i<n;i++)
        {
            upd_u[i] = u[i]-(c*dt/dx)*(u[i]-u[i-1]);
        }

        
    }
    return upd_u;
}

// QUICK scheme
double* QUICK_scheme(int n,double* &u,double dx, double dt, double c,int proc_id,int num_procs, double *fwd_pts, double* bwd_pts)
{
    double* upd_u = new double[n];

    if(proc_id ==0)
    {
        upd_u[1] = u[1]-(c*dt/dx)*(u[1]-u[0]);

        for(int i=2;i<n-1;i++)
        {
            upd_u[i] = u[i]-(c*dt/dx)*((3.0/8.0)*u[i]-(7.0/8.0)*u[i-1]+(1.0/8.0)*u[i-2]+(3.0/8.0)*u[i+1]);
        }
        upd_u[n-1] = u[n-1]-(c*dt/dx)*((3.0/8.0)*u[n-1]-(7.0/8.0)*u[n-2]+(1.0/8.0)*u[n-3]+(3.0/8.0)*fwd_pts[0]);
    }
    else
    {
        double prev_proc_last_pt = bwd_pts[(sizeof(bwd_pts))/sizeof(bwd_pts[0])];

        double prev_proc_scnd_last_pt = bwd_pts[(sizeof(bwd_pts))/sizeof(bwd_pts[0])-1];

        upd_u[0] = u[0]- (c*dt/dx)*((3.0/8.0)*u[0]-(7.0/8.0)*prev_proc_last_pt+(1.0/8.0)*prev_proc_scnd_last_pt+(3.0/8.0)*u[1]);

        upd_u[1] = u[1]- (c*dt/dx)*((3.0/8.0)*u[1]-(7.0/8.0)*u[0]+(1.0/8.0)*prev_proc_last_pt+(3.0/8.0)*u[2]);
        for(int i=2;i<n-1;i++)
        {
            
            upd_u[i] = u[i]-(c*dt/dx)*((3.0/8.0)*u[i]-(7.0/8.0)*u[i-1]+(1.0/8.0)*u[i-2]+(3.0/8.0)*u[i+1]);
            
        }

        if(proc_id!=num_procs-1)
        {
            upd_u[n-1] = u[n-1]-(c*dt/dx)*((3.0/8.0)*u[n-1]-(7.0/8.0)*u[n-2]+(1.0/8.0)*u[n-3]+(3.0/8.0)*fwd_pts[0]);
        }
        else
        {
            upd_u[n-1] = 0;
        }
    }
    return upd_u;
    
    
}

// Function to write the vector of double* values as a csv file. (Also adds a column of the time step )

void wave_data_to_csv(int n, double** data, double dt, double dx, double left_pt, double right_pt, double max_time, const std::string& filename)
{
    const string path = "Solution/Q1b/"+filename;

    std::ofstream outputFile(path);

    if (!outputFile.is_open()) 
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    
    //Header row with name and the x positions
    outputFile << "Time";
    for (double i = left_pt; i <= right_pt; i+=dx) 
    {
        outputFile << "," << i ;
    }
    outputFile << "\n";

    

    int numTimeSteps = round(max_time / dt);      //Total number of time steps


    for (int step = 0; step <= numTimeSteps; step++) 
    {
        double currentTime = step * dt;
        if(currentTime == 0 || currentTime== 0.5  || currentTime == 1.0)
        {
            outputFile << currentTime;                          //Adding the time step currentTime to first column
            for (int col = 0; col < n; col++) 
            {
            
                outputFile << "," << (data[step][col]);         //Adding the wave data at time step currentTime to a row
            }
            outputFile << std::endl;
        }
        
    }
    outputFile.close();
    
    cout<<"Successfully wrote the data to a csv file !"<<endl<<endl;
}


int main(int argc, char** argv)
{
    //define variables for MPI communication
    int i, j, proc_id, num_pros, tag=100;
    MPI_Status status;
    

    //Initialize MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_pros);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    //Define all the variables
    double L = 2;
    double left_pt = 0;
    double right_pt = L;
    double mid_pt = 0.5;
    double dx = 0.002;
    double dt = 0.0001;

    int n = (right_pt - left_pt)/dx;    

    double c = 1.0;



    double local_start = ((right_pt-left_pt)/num_pros) + (proc_id-1)*((right_pt-left_pt)/num_pros);
    int local_n = n/num_pros;


    int shdw_pt_fwd = 1;
    int shdw_pt_bwd = 2;
    double* recv_bwd_shdw_pts = new double[shdw_pt_bwd];
    double* recv_fwd_shdw_pts = new double[shdw_pt_fwd];

    

    
    
    double max_time= 1;



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                              Analytical solution
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(proc_id==0)
    {
        double* u0 = initialize_wave(n,left_pt,mid_pt,right_pt,dx);
        double* local_u0 = new double[local_n];
        cout<<"Performing the analytical solution calculation"<<endl;

        double** anal_solution_matrix = new double*[int(max_time/dt)];
        for(int i=0; i<=int(max_time/dt);i++)
        {
            anal_solution_matrix[i] = new double[n]{};
        }
        
        //Define the solution variables
        double* anal_sol = u0;

        double* anal_recv_arr = new double[local_n]{};
        
        //Add the initial wave data to solution matrices
        
        for(int i=0;i<n;i++)
        {
            anal_solution_matrix[0][i] = u0[i];
        }
        

        for(int i=0;i<local_n;i++)
        {
            local_u0[i] = u0[i];
        }

        //sending the array to each of the processor
        for(int i=1;i<num_pros;i++)
        {
            MPI_Send( &u0[local_n*i], local_n , MPI_DOUBLE ,i, tag, MPI_COMM_WORLD);
        }


        for(double anal_time = dt; anal_time<=max_time;anal_time+=dt)
        {
            
            local_u0 = analytical_solution(local_n,local_start,left_pt,mid_pt,right_pt,dx,c,anal_time);
            
            for(int i=0;i<local_n;i++)
            {
                anal_sol[i] = local_u0[i];
            }
    
            for(int i=1;i<num_pros;i++)
            {
                
                MPI_Recv( anal_recv_arr, local_n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                
                for(int j=0;j<local_n;j++)
                {
                    anal_sol[i*local_n+j] = anal_recv_arr[j];

                }

            }

            for(int i=0;i<n;i++)
            {
                anal_solution_matrix[int(round(anal_time/dt))][i] = anal_sol[i];
            }
        }

        cout<<"Successfully completed the analytical solution calculation"<<endl;
        
        cout<<"Writing the analytical solution data as a csv file..."<<endl;
        string file_name = "1D_Wave_MPI_"+to_string(num_pros)+"_analytical.csv";
        wave_data_to_csv(n,anal_solution_matrix,dt,dx,left_pt,right_pt,max_time,file_name);

        delete[] anal_recv_arr;
        for(int i=0;i<round(max_time/dt);i++)
        {
            delete[] anal_solution_matrix[i];
        }
        delete[] anal_solution_matrix;
        delete[] local_u0;
        delete[] u0;
    }

    if(proc_id!=0)
    {   
        double* local_u0 = new double[local_n];
        MPI_Recv( local_u0 ,local_n , MPI_DOUBLE , 0 , tag , MPI_COMM_WORLD , &status);

        for(double anal_time=dt;anal_time<=max_time;anal_time+=dt)
        {

            local_u0 = analytical_solution(local_n,local_start,left_pt,mid_pt,right_pt,dx,c,anal_time);

            MPI_Send( local_u0,local_n,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

        }
        delete[] local_u0;
    }




    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                              UPWIND SCHEME
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////






    if(proc_id==0)
    {
        double* u0 = initialize_wave(n,left_pt,mid_pt,right_pt,dx);
        double* local_u0 = new double[local_n];
        cout<<"Performing the upwind scheme solution calculation"<<endl;

        double** upwind_solution_matrix = new double*[int(max_time/dt)];
        for(int i=0; i<=int(max_time/dt);i++)
        {
            upwind_solution_matrix[i] = new double[n]{};
        }
        
        //Define the solution variables
        double* upwind_sol = u0;

        double* upwind_recv_arr = new double[local_n]{};
        
        //Add the initial wave data to solution matrices
        
        for(int i=0;i<n;i++)
        {
            upwind_solution_matrix[0][i] = u0[i];
        }
        

        for(int i=0;i<local_n;i++)
        {
            local_u0[i] = u0[i];
        }

        //sending the array to each of the processor
        for(int i=1;i<num_pros;i++)
        {
            MPI_Send( &u0[local_n*i], local_n , MPI_DOUBLE ,i, tag, MPI_COMM_WORLD);
        }


        for(double upwind_time = dt; upwind_time<=max_time;upwind_time+=dt)
        {
            
            //Sending shadow points
            MPI_Send( &local_u0[local_n-shdw_pt_bwd], shdw_pt_bwd , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);                  //Last two points to succeeding processor
            
            local_u0 = upwind_scheme(local_n,local_u0,dx,dt,c,proc_id,recv_bwd_shdw_pts);
            
            for(int i=0;i<local_n;i++)
            {
                upwind_sol[i] = local_u0[i];
            }
    
            for(int i=1;i<num_pros;i++)
            {
                
                MPI_Recv( upwind_recv_arr, local_n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                
                for(int j=0;j<local_n;j++)
                {
                    upwind_sol[i*local_n+j] = upwind_recv_arr[j];
                }

            }

            


            for(int i=0;i<n;i++)
            {
                upwind_solution_matrix[int(round(upwind_time/dt))][i] = upwind_sol[i];
            }
        }
        
        cout<<"Successfully completed the upwind scheme solution calculation"<<endl;

        cout<<"Writing the upwind solution data as a csv file..."<<endl;
        string file_name = "1D_Wave_MPI_"+to_string(num_pros)+"_upwind.csv";
        wave_data_to_csv(n,upwind_solution_matrix,dt,dx,left_pt,right_pt,max_time,file_name);

        delete[] upwind_recv_arr;
        for(int i=0;i<round(max_time/dt);i++)
        {
            delete[] upwind_solution_matrix[i];
        }
        delete[] upwind_solution_matrix;
        delete[] local_u0;
        delete[] u0;
    }

    if(proc_id!=0)
    {   
        
        double* local_u0 = new double[local_n];
        MPI_Recv( local_u0 ,local_n , MPI_DOUBLE , 0 , tag , MPI_COMM_WORLD , &status);

        for(double upwind_time=dt;upwind_time<=max_time;upwind_time+=dt)
        {

            if(proc_id!=num_pros-1)
            {   

                //Sending shadow points
                MPI_Send( &local_u0[local_n-shdw_pt_bwd], shdw_pt_bwd , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);            //Last two points to succeeding processor

                //Receiving shadow points
                MPI_Recv( recv_bwd_shdw_pts ,shdw_pt_bwd , MPI_DOUBLE , proc_id-1 , tag , MPI_COMM_WORLD , &status);            //Last two points from preceding processor
                
            }

            else
            {
                //Receiving shadow points
                MPI_Recv( recv_bwd_shdw_pts ,shdw_pt_bwd , MPI_DOUBLE , proc_id-1 , tag , MPI_COMM_WORLD , &status);            //Last two points from preceding processor

            }

            local_u0 = upwind_scheme(local_n,local_u0,dx,dt,c,proc_id,recv_bwd_shdw_pts);

            MPI_Send( local_u0,local_n,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
        }
        delete[] local_u0;
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                              QUICK SCHEME
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(proc_id==0)
    {
        double* u0 = initialize_wave(n,left_pt,mid_pt,right_pt,dx);
        double* local_u0 = new double[local_n];
        cout<<"Performing the QUICK scheme solution calculation"<<endl;
        double** QUICK_solution_matrix = new double*[int(max_time/dt)];
        for(int i=0; i<=int(max_time/dt);i++)
        {
            QUICK_solution_matrix[i] = new double[n]{};
        }

        
        //Define the solution variables
        double* QUICK_sol = u0;

        double* QUICK_recv_arr = new double[local_n]{};
        
        //Add the initial wave data to solution matrices
        for(int i=0;i<n;i++)
        {
            QUICK_solution_matrix[0][i] = u0[i];
        }
        

        for(int i=0;i<local_n;i++)
        {
            local_u0[i] = u0[i];
        }

        //sending the array to each of the processor
        for(int i=1;i<num_pros;i++)
        {
            MPI_Send( &u0[local_n*i], local_n , MPI_DOUBLE ,i, tag, MPI_COMM_WORLD);
        }


        for(double QUICK_time = dt; QUICK_time<=max_time;QUICK_time+=dt)
        {
            
            //Sending shadow points
            
            MPI_Send( &local_u0[local_n-shdw_pt_bwd], shdw_pt_bwd , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);                  //Last two points to succeeding processor
            
            //Receive shadow points
            
            MPI_Recv( recv_fwd_shdw_pts ,shdw_pt_fwd , MPI_DOUBLE , proc_id+1 , tag , MPI_COMM_WORLD , &status);            //First point from succeeding processor
        
            
            local_u0 = QUICK_scheme(local_n,local_u0,dx,dt,c,proc_id,num_pros,recv_fwd_shdw_pts,recv_bwd_shdw_pts);

            for(int i=0;i<local_n;i++)
            {
                QUICK_sol[i] = local_u0[i];
            }
    
            for(int i=1;i<num_pros;i++)
            {
                
                MPI_Recv( QUICK_recv_arr, local_n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                
                for(int j=0;j<local_n;j++)
                {
                    QUICK_sol[i*local_n+j] = QUICK_recv_arr[j];
                }

            }


            for(int i=0;i<n;i++)
            {
                QUICK_solution_matrix[int(round(QUICK_time/dt))][i] = QUICK_sol[i];
            }

        }
        cout<<"Successfully completed the QUICK scheme solution calculation"<<endl;

        cout<<"Writing the QUICK solution data as a csv file..."<<endl;
        string file_name = "1D_Wave_MPI_"+to_string(num_pros)+"_QUICK.csv";
        wave_data_to_csv(n,QUICK_solution_matrix,dt,dx,left_pt,right_pt,max_time,file_name);

        delete[] QUICK_recv_arr;
        for(int i=0;i<round(max_time/dt);i++)
        {
            delete[] QUICK_solution_matrix[i];
        }
        delete[] QUICK_solution_matrix;
        delete[] local_u0;
        delete[] u0;
    }

    if(proc_id!=0)
    {   
        
        double* local_u0 = new double[local_n];
        MPI_Recv( local_u0 ,local_n , MPI_DOUBLE , 0 , tag , MPI_COMM_WORLD , &status);

        for(double QUICK_time=dt;QUICK_time<=max_time;QUICK_time+=dt)
        {

            if(proc_id!=num_pros-1)
            {   

                //Sending shadow points
                MPI_Send( &local_u0[local_n-shdw_pt_bwd], shdw_pt_bwd , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);            //Last two points to succeeding processor

                //Receiving shadow points
                MPI_Recv( recv_fwd_shdw_pts ,shdw_pt_fwd , MPI_DOUBLE , proc_id+1 , tag , MPI_COMM_WORLD , &status);            //First point from succeeding processor

                //Sending shadow points
                MPI_Send( local_u0, shdw_pt_fwd , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD);                                 //First one point to preceding processor
                
                //Receiving shadow points
                MPI_Recv( recv_bwd_shdw_pts ,shdw_pt_bwd , MPI_DOUBLE , proc_id-1 , tag , MPI_COMM_WORLD , &status);            //Last two points from preceding processor
                
            }

            else
            {
                //Sending shadow points
                MPI_Send( local_u0, shdw_pt_fwd , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD);                                  //First one point to preceding processor

                //Receiving shadow points
                MPI_Recv( recv_bwd_shdw_pts ,shdw_pt_bwd , MPI_DOUBLE , proc_id-1 , tag , MPI_COMM_WORLD , &status);            //Last two points from preceding processor

            }

            local_u0 = QUICK_scheme(local_n,local_u0,dx,dt,c,proc_id,num_pros,recv_fwd_shdw_pts,recv_bwd_shdw_pts);

            MPI_Send( local_u0,local_n,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

        }
        delete[] local_u0;
    }




    
    delete[] recv_fwd_shdw_pts; // Free memory for shadow point arrays
    delete[] recv_bwd_shdw_pts;

    MPI_Finalize();
}