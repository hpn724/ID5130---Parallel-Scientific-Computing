#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

int main() 
{
    // Compile MPI program
    string compile_commandQ1 = "mpiCC -o Tut7_Q1 Tut7_Q1.cpp";
    string compile_commandQ2 = "mpiCC -o Tut7_Q2 Tut7_Q2.cpp";
    string compile_commandQ3 = "mpiCC -o Tut7_Q3 Tut7_Q3.cpp";
    
    
    int compile_status_1 = system(compile_commandQ1.c_str());

    if (compile_status_1 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }

    int compile_status_2 = system(compile_commandQ2.c_str());

    if (compile_status_2 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }

    int compile_status_3 = system(compile_commandQ3.c_str());

    if (compile_status_3 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }

    // Run MPI program with desired number of processors

    string num_procs_str = "4"; // Number of processors as a string
    
    string run_command_1 = "mpiexec -np " + num_procs_str + " ./Tut7_Q1";

    int run_status_1 = system(run_command_1.c_str());

    if (run_status_1 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }

    string run_command_2 = "mpiexec -np " + num_procs_str + " ./Tut7_Q2";

    int run_status_2 = system(run_command_2.c_str());

    if (run_status_2 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }

    string run_command_3 = "mpiexec -np " + num_procs_str + " ./Tut7_Q3";

    int run_status_3 = system(run_command_3.c_str());

    if (run_status_3 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }

    return 0;
}
