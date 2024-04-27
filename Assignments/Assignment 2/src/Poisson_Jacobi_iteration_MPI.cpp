#include<iostream>
#include<cmath>
#include<fstream>
#include<string.h>
#include "mpi.h"
#include <sstream>
#include <iomanip>
#include<vector>

using namespace std;

#define pi 3.14159265358

void show_matrix(int m, int n, double** mat)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }
}

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

double q_func(double x, double y)
{
    return pow(x,2)+pow(y,2);
}

void initialize_grid(double delta,int n_x,int n_y, double x_start, double y_start, double* &x, double* &y)
{
    
    for(int i=0;i<n_x;i++)
    {   
        x[i] = x_start + i*delta;
    }

    for(int i=0;i<n_y;i++)
    {   
        y[i] = y_start + i*delta;
    }
}

void initialize_solution(int n_x, int n_y, double *x, double *y, double** &phi, double** &phi_upd, double** &q,int proc_id)
{
    phi = new double*[n_x];
    phi_upd = new double*[n_x];
    q = new double*[n_x];

    for(int i=0;i<n_x;i++)
    {
        phi[i] = new double[n_y];
        phi_upd[i] = new double[n_y];
        q[i] = new double[n_y];
        for(int j=0; j<n_y;j++)
        {
            q[i][j] = q_func(x[i],y[j]);
            if(i==0 && proc_id==0)
            {
                phi[i][j] = sin(2*pi*y[j]);
                phi_upd[i][j] = sin(2*pi*y[j]);
            }
            else
            {
                phi[i][j] = 0;
                phi_upd[i][j] = 0;
            }
        }
    }
}

double** jacobi_iteration(double delta,int local_n_x, int n_y,double *y, double** local_phi_upd, double** q, int num_pros, int proc_id, double* shdw_pts_up, double* shdw_pts_dwn)
{
    
    for(int j=1;j<n_y-1;j++)
    {
        if(proc_id==0)
        {
            local_phi_upd[0][j] = sin(2*pi*y[j]);
        }
        else
        {
            local_phi_upd[0][j] = 0.25*(local_phi_upd[1][j]+shdw_pts_up[j]+local_phi_upd[0][j+1]+local_phi_upd[0][j-1] + pow(delta,2)*q[0][j]);
        }
    }

    for(int j=1;j<n_y-1;j++)
    {
        for(int i=1;i<local_n_x-1;i++)
        {

            local_phi_upd[i][j] = 0.25*(local_phi_upd[i+1][j]+local_phi_upd[i-1][j]+local_phi_upd[i][j+1]+local_phi_upd[i][j-1] + pow(delta,2)*q[i][j]);
            
        }
    }



    for(int j=1;j<n_y-1;j++)
    {
        if(proc_id==num_pros-1)
        {
            local_phi_upd[local_n_x-1][j] = (4*local_phi_upd[local_n_x-2][j] -local_phi_upd[local_n_x-3][j])/3;
        }
        else
        {
            local_phi_upd[local_n_x-1][j] = 0.25*(shdw_pts_dwn[j]+local_phi_upd[local_n_x-2][j]+local_phi_upd[local_n_x-1][j+1]+local_phi_upd[local_n_x-1][j-1] + pow(delta,2)*q[local_n_x-1][j]);
        }
        
    }

    return local_phi_upd;

}

double norm(int n_x, int n_y, double **phi, double** phi_upd)
{
    double err=0;
    for(int i=0;i<n_x;i++)
    {
        for(int j=0;j<n_y;j++)
        {
            err += pow(phi_upd[i][j]-phi[i][j],2);
        }
    }
    return sqrt(err);
}

void write_to_csv(double delta,int n_x, int n_y, double *x, double *y, double** phi_upd, string filename)
{
    
    string path;
    if(delta == 0.005)
    {
        path = "Solution/Q2d/"+filename;
    }
    else
    {
        path = "Solution/Q2b/"+filename;
    }

    std::ofstream outputFile(path);

    if (!outputFile.is_open()) 
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    
    
    outputFile << "Coordinates";
    for (int i = 0; i < n_y; i++) 
    {
        outputFile << "," << y[i] ;
    }
    outputFile << "\n";

    
    for (int i=0;i<n_x;i++) 
    {
        outputFile << x[i];
        for (int j = 0; j < n_y; j++) 
        {
            outputFile << "," << (phi_upd[i][j]);  
        }
        outputFile << std::endl;
        
    }

    outputFile.close();
    cout<<"Successfully wrote the data to a csv file !"<<endl<<endl;
}

void vector_csv(double delta, vector<long int> iter, vector<double> errors, string filename)
{
    string path;
    if(delta == 0.005)
    {
        path = "Solution/Q2d/"+filename;
    }
    else
    {
        path = "Solution/Q2b/"+filename;
    }
    

    std::ofstream outputFile(path);

    if (!outputFile.is_open()) 
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    
    
    outputFile << "Iteration,Error"<<endl;
    for (int i = 0; i < errors.size(); i++) 
    {
        outputFile <<iter[i]<< "," << errors[i]<<endl ;
    }


    outputFile.close();
    cout<<"Successfully wrote iteration and error data to a csv file !"<<endl<<endl;
}

void time_data(double time_taken, double serial_time, int num_procs, string filename)
{
    const string filepath = "Solution/Q2d/"+filename;

    std::ofstream outputFile(filepath, std::ios_base::app);
    
    if (!outputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    double f = serial_time/time_taken;

    double psi = 1/(f+((1-f)/num_procs));
    
    
    outputFile << num_procs<<","<<psi << std::endl;
    
    outputFile.close();
}



int main(int argc, char** argv)
{
    int proc_id, num_pros, tag =100;
    int tag_y = 101;
    int tag_phi = 1002;
    int tag_q = 0;
    
    MPI_Status status;

    //Initialize MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_pros);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);


    double time_taken = MPI_Wtime();
    
    double delta = 0.005;

    double x_start = -1;
    double y_start = -1;
    double x_end = 1;
    double y_end = 1;

    int n_x = round((x_end-x_start)/delta)+1;
    int n_y = round((y_end-y_start)/delta)+1;

    int local_n_x = n_x/num_pros;
    
    if(proc_id==0 && n_x%local_n_x!=0)
    {
        local_n_x+=1;
    }

    double *x = new double[n_x];
    double *y = new double[n_y];

    initialize_grid(delta,n_x,n_y,x_start,y_start,x,y);

    double *local_x = new double[local_n_x];
    double *local_y = new double[local_n_x];
    for(int i=0;i<local_n_x;i++)
    {
        if(proc_id==0)
        {
            local_x[i] = x[i];
            local_y[i] = y[i];
        }
        else
        {
            local_x[i] = x[proc_id*local_n_x+1+i];
            local_y[i] = y[proc_id*local_n_x+1+i];
        }
    }

    double **local_phi = nullptr;
    double **local_phi_upd = nullptr;
    double **local_q = nullptr;

    initialize_solution(local_n_x, n_y, local_x, y, local_phi, local_phi_upd, local_q,proc_id);
    

    double err_thresh = 0.0001;

    int iter_lim = pow(10,7);

    long int iter_count=0;

    double curr_err = 100;

    double* shdw_pts_up = new double[n_y]{};
    double* shdw_pts_dwn = new double[n_y]{};

    double send_arr[local_n_x][n_y];

    if(proc_id==0)
    {
        double **phi = nullptr;
        double **phi_upd = nullptr;
        double **q = nullptr;

        double recv_arr[local_n_x-1][n_y];


        vector<double> errors;
        vector<long int> iterations;
        

        initialize_solution(n_x, n_y, x, y, phi, phi_upd, q, proc_id);

        double serial_time = 0;
        
        double temp_time ;

        while(abs(curr_err) > err_thresh && iter_count < iter_lim)
        {
            
            MPI_Send( &(local_phi_upd[local_n_x-1][0]), n_y , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);

            MPI_Recv( &(shdw_pts_dwn[0]), n_y , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD, &status);

            iter_count += 1;
            for(int i=1;i<num_pros;i++)
            {
                MPI_Send( &iter_count, 1 , MPI_LONG ,i, tag, MPI_COMM_WORLD);
            }

            local_phi_upd = jacobi_iteration(delta, local_n_x, n_y, y ,local_phi_upd,local_q, num_pros, proc_id, shdw_pts_dwn, shdw_pts_up);



            for(int i=0;i<local_n_x;i++)
            {
                for(int j=0;j<n_y;j++)
                {
                    phi_upd[i][j] = local_phi_upd[i][j];
                }
            }
            
            for(int i=1;i<num_pros;i++)
            {
                MPI_Recv(&(recv_arr[0][0]), (local_n_x-1)*n_y, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                
                for(int j=0;j<local_n_x-1;j++)
                {
                    for(int k=0;k<n_y;k++)
                    {
                        
                        phi_upd[i*local_n_x-(i-1)+j][k] = recv_arr[j][k];
                        
                    }
                }

            }
    
            temp_time = MPI_Wtime();
            curr_err = norm(n_x,n_y,phi,phi_upd);

            

            for(int i=1;i<num_pros;i++)
            {
                MPI_Send(&curr_err, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            }

            for(int i=0;i<n_x;i++)
            {
                for(int j=0;j<n_y;j++)
                {
                    phi[i][j]=phi_upd[i][j];
                }
            }
            cout<< "Error : "<<curr_err<<" \t \t  iteration : "<<iter_count<<endl;

            iterations.push_back(iter_count);
            errors.push_back(curr_err);

            temp_time = MPI_Wtime() - temp_time;

            serial_time += temp_time;

        }


        if(curr_err<err_thresh)
        {
            cout<<"Iteration completed. Solution has converged !"<<endl;
        }
        else
        {
            cout<<"Iteration completed. Solution failed to converge"<<endl;
        }

        time_taken = MPI_Wtime() - time_taken;

        
        cout<<"Time taken for the program to run : "<<time_taken<<endl;


        std::stringstream ss;
        
        ss << std::fixed << std::setprecision(3) << delta;
        
        std::string delta_s = ss.str();

        cout<<"Writing the final matrix as a csv file"<<endl;
        string filename_mat = "Poisson_Jacobi_MPI_"+ to_string(num_pros) +"_"+delta_s+".csv";
        write_to_csv(delta,n_x,n_y,x,y,phi_upd,filename_mat);

        string filename_vec = "Poisson_Jacobi_error_MPI_"+ to_string(num_pros) +"_"+delta_s+".csv";
        vector_csv(delta,iterations,errors,filename_vec);

        if(delta == 0.005)
        {
            string filename_time = "Poisson_Jacobi_speedup.csv";
            time_data(time_taken,serial_time,num_pros,filename_time);
        }


        for(int i=0;i<n_x;i++)
        {
            delete[] phi[i];
            delete[] phi_upd[i];
            delete[] q[i];

        }
        delete[] phi;
        delete[] phi_upd;
        delete[] q;

    }

    else if(proc_id!=0)
    {
        
        
        while(abs(curr_err) > err_thresh && iter_count < iter_lim)
        {
            
            if(proc_id!=num_pros-1)
            {
                MPI_Send( &(local_phi_upd[local_n_x-1][0]), n_y , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD);

                MPI_Recv( &(shdw_pts_dwn[0]), n_y , MPI_DOUBLE ,proc_id+1, tag, MPI_COMM_WORLD, &status);
                
                MPI_Send( &(local_phi_upd[0][0]), n_y , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD);

                MPI_Recv( &(shdw_pts_up[0]), n_y , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD, &status);
                
            }

            else
            {
                MPI_Send( &(local_phi_upd[0][0]), n_y , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD);

                MPI_Recv( &(shdw_pts_up[0]), n_y , MPI_DOUBLE ,proc_id-1, tag, MPI_COMM_WORLD, &status);
            

            }
            
            MPI_Recv( &iter_count, 1,MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
            
            local_phi_upd = jacobi_iteration(delta, local_n_x, n_y, y ,local_phi_upd,local_q, num_pros, proc_id, shdw_pts_dwn, shdw_pts_up);
            

            for(int i=0;i<local_n_x;i++)
            {
                for(int j=0;j<n_y;j++)
                {
                    send_arr[i][j] = local_phi_upd[i][j];
                }
            }
            
            
            MPI_Send( &(send_arr[0][0]), local_n_x*n_y, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );

            MPI_Recv( &curr_err, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

        }
        
        
    }
    
        
    delete[] x;
    delete[] y;

    for(int i=0;i<local_n_x;i++)
    {
        delete[] local_phi[i];
        delete[] local_phi_upd[i];
        delete[] local_q[i];
    }


    delete[] local_phi;
    delete[] local_phi_upd;
    delete[] local_q;

    delete[] local_y;
    delete[] local_x;

    MPI_Finalize();
    
    return 0;


}