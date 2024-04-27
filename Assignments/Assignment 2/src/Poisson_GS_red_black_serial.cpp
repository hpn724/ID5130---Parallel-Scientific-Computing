#include<iostream>
#include<cmath>
#include<fstream>
#include<string.h>
#include<sstream>
#include<iomanip>

using namespace std;

#define pi 3.14159265358

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

void initialize_solution(int n_x, int n_y, double *x, double *y, double** &phi, double** &phi_upd, double** &q)
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
            if(i==0)
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

double** Gauss_Seidel_red_black(double delta,int n_x, int n_y,double *y, double** phi_upd, double** q)
{
    for(int j=0;j<n_y;j++)
    {
        phi_upd[0][j] = sin(2*pi*y[j]);
    }

    for(int j=1;j<n_y-1;j++)
    {
        for(int i=1;i<n_x-1;i++)
        {
            if((i+j)%2==0)
            {
                phi_upd[i][j] = 0.25*(phi_upd[i+1][j]+phi_upd[i-1][j]+phi_upd[i][j+1]+phi_upd[i][j-1] + pow(delta,2)*q[i][j]);
            }
        }
    }

    for(int j=1;j<n_y-1;j++)
    {
        for(int i=1;i<n_x-1;i++)
        {
            if((i+j)%2!=0)
            {
                phi_upd[i][j] = 0.25*(phi_upd[i+1][j]+phi_upd[i-1][j]+phi_upd[i][j+1]+phi_upd[i][j-1] + pow(delta,2)*q[i][j]);
            }
        }
    }

    for(int j=0;j<n_y;j++)
    {
        phi_upd[n_x-1][j] = (4*phi_upd[n_x-2][j] -phi_upd[n_x-3][j])/3;
    }

    return phi_upd;

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

void write_to_csv(double delta, int n_x, int n_y, double *x, double *y, double** phi_upd, string filename)
{
        const string path = "../Solution/Q2c/"+filename;

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


int main()
{
    double delta = 0.1;

    double x_start = -1;
    double y_start = -1;
    double x_end = 1;
    double y_end = 1;

    int n_x = round((x_end-x_start)/delta)+1;
    int n_y = round((y_end-y_start)/delta)+1;

    double *x = new double[n_x];
    double *y = new double[n_y];

    initialize_grid(delta,n_x,n_y,x_start,y_start,x,y);

    double **phi, **phi_upd, **q;

    initialize_solution(n_x,n_y,x,y,phi,phi_upd,q);

    double err_thresh = 0.0001;

    int iter_lim = pow(10,7);

    long int iter_count=0;

    double curr_err = 100;

    while(abs(curr_err) > err_thresh && iter_count < iter_lim)
    {
        iter_count += 1;
        phi_upd = Gauss_Seidel_red_black(delta,n_x,n_y,y,phi_upd,q);
        curr_err = norm(n_x,n_y,phi,phi_upd);

        for(int i=0;i<n_x;i++)
        {
            for(int j=0;j<n_y;j++)
            {
                phi[i][j]=phi_upd[i][j];
            }
        }
        cout<< "Error : "<<curr_err<<" \t \t  iteration : "<<iter_count<<endl;

    }


    if(curr_err<err_thresh)
    {
        cout<<"Iteration completed. Solution has converged !";
    }
    else
    {
        cout<<"Iteration completed. Solution failed to converge";
    }

    std::stringstream ss;
        
    ss << std::fixed << std::setprecision(3) << delta;
    
    std::string delta_s = ss.str();

    cout<<"Writing the final matrix as a csv file"<<endl;

    string file_name = "Poisson_GS_red_black_serial" + delta_s+".csv";

    write_to_csv(delta,n_x,n_y,x,y,phi_upd,file_name);
    


    delete[] x;
    delete[] y;
    
    for(int i=0;i<n_x;i++)
    {
        delete[] phi[i];
        delete[] phi_upd[i];
        delete[] q[i];
    }

    delete[] phi;
    delete[] phi_upd;
    delete[] q;

    return 0;

}