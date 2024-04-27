#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

using namespace std;

#define pi 3.14159265358

// function to show the vector (to verify the output if the o/p is correct or not)
void show_vector(int n, double* x)
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
double* initialize_wave(int n,double left_pt, double mid_pt, double right_pt, double dx)
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
double* analytical_solution(int n,double left_pt, double mid_pt, double right_pt, double dx, double c, double t)
{
    double* u = new double[n];
    double curr_x = left_pt;
    for(int i=0;i<n;i++)
    {
        u[i] = initial_wave_function((curr_x-(c*t)),left_pt,mid_pt,right_pt);
        curr_x+=dx;
    }
    return u;
}

// Upwind scheme
double* upwind_scheme(int n,double* &u, double dx, double dt, double c)
{
    double* upd_u = new double[n];
    for(int i=1;i<n;i++)
    {

        upd_u[i] = u[i]-(c*dt/dx)*(u[i]-u[i-1]);
        
    }
    return upd_u;
}

// QUICK scheme
double* QUICK_scheme(int n,double* &u,double dx, double dt, double c)
{
    double* upd_u = new double[n];
    for(int i=2;i<n;i++)
    {
        upd_u[i] = u[i]-(c*dt/dx)*((3.0/8.0)*u[i]-(7.0/8.0)*u[i-1]+(1.0/8.0)*u[i-2]+(3.0/8.0)*u[i+1]);
    }
    upd_u[1] = u[1]-(c*dt/dx)*(u[1]-u[0]);
    return upd_u;
}

// Function to write the vector of double* values as a csv file. (Also adds a column of the time step )
void wave_data_to_csv(int n, const std::vector<double*>& data, double dt, double dx, double left_pt, double right_pt, double max_time, const std::string& filename)
{
    const string path = "../Solution/Q1a/"+filename;

    std::ofstream outputFile(path);

    if (!outputFile.is_open()) 
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    
    //Header row with name and the x positions
    outputFile << "Time";
    for (double i = left_pt; i < right_pt; i+=dx) 
    {
        outputFile << "," << i ;
    }
    outputFile << "\n";

    

    int numTimeSteps = static_cast<int>(max_time / dt);      //Total number of time steps


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


int main()
{

    //Define all the variables
    double L = 2;
    double left_pt = 0;
    double right_pt = L;
    double mid_pt = 0.5;
    double dx = 0.002;
    double dt = 0.0001;

    int n = (right_pt - left_pt)/dx;    

    double c = 1.0;


    //Define the solution variables
    double* u0 = initialize_wave(n,left_pt,mid_pt,right_pt,dx);
    double* upwind_sol = u0;
    double* QUICK_sol = u0;
    double* anal_sol = u0;

    double time_upwind = 0.0;
    double time_QUICK = 0.0;
    double max_time= 1.0;

    vector<double*> upwind_solution_matrix;
    vector<double*> QUICK_solution_matrix;
    vector<double*> analytical_sol_matrix;

    //Add the initial wave data to solution matrices
    upwind_solution_matrix.push_back(upwind_sol);
    QUICK_solution_matrix.push_back(QUICK_sol);
    analytical_sol_matrix.push_back(anal_sol);

    

    //Perform calculation for each time step
    cout<<"Performing upwind scheme, QUICK scheme and analytical solution calculations..."<<endl;

    for(double time=0;time<=max_time;time+=dt)
    {
        if(time<max_time)
        {
            upwind_sol = upwind_scheme(n, upwind_sol, dx, dt, c);
            upwind_solution_matrix.push_back(upwind_sol);

            QUICK_sol = QUICK_scheme(n, QUICK_sol, dx, dt, c);
            QUICK_solution_matrix.push_back(QUICK_sol);

            anal_sol = analytical_solution( n, left_pt,  mid_pt, right_pt,  dx,  c,  time);
            analytical_sol_matrix.push_back(anal_sol);
        }
    }

    cout<<"Successfully completed calculating the values and filling up of the solution matrices !!"<<endl;

    //Write the data to csv files
    cout<<endl<<"Writing upwind scheme data to a csv file"<<endl;
    wave_data_to_csv(n,upwind_solution_matrix,dt,dx,left_pt,right_pt,max_time,"1D_Wave_serial_upwind.csv");
    cout<<"Writing QUICK scheme data to a csv file"<<endl;
    wave_data_to_csv(n,QUICK_solution_matrix,dt,dx,left_pt,right_pt,max_time,"1D_Wave_serial_QUICK.csv");
    cout<<"Writing analytical solution data to a csv file"<<endl;
    wave_data_to_csv(n,analytical_sol_matrix,dt,dx,left_pt,right_pt,max_time,"1D_Wave_serial_analytical.csv");
    

    //Free up the used memory
    delete[] u0;
    delete[] upwind_sol;
    delete[] QUICK_sol;
    delete[] anal_sol;

    // Clear and free memory for vectors of solution matrices
    for (auto ptr : upwind_solution_matrix)
    {
        delete[] ptr;
    }
    upwind_solution_matrix.clear();

    for (auto ptr : QUICK_solution_matrix)
    {
        delete[] ptr;
    }
    QUICK_solution_matrix.clear();

    for (auto ptr : analytical_sol_matrix)
    {
        delete[] ptr;
    }
    analytical_sol_matrix.clear();

    return 0;

}