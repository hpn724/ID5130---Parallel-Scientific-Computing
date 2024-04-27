#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

int main() 
{
    /*

    //////////////////////////////////////////////////////////
    //                      Question 1                      //
    //////////////////////////////////////////////////////////

    string compile_commandQ1 = "mpiCC -o Tut8_Q1 Tut8_Q1.cpp";
    
    int compile_status_1 = system(compile_commandQ1.c_str());

    if (compile_status_1 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }
    
    vector<string> num_procs = {"2","4","6"};

    for(auto procs : num_procs)
    {
        string run_command_1 = "mpiexec -np " + procs + " ./Tut8_Q1";
        cout<<"**********************************************"<<endl;
        cout<<"Running Question 1 with "<<procs<<" number of processors"<<endl;
        cout<<"**********************************************"<<endl;
        int run_status_1 = system(run_command_1.c_str());

        if (run_status_1 != 0) 
        {
            cerr << "Error running MPI program." << endl;
            return 1;
        }
    }


    //////////////////////////////////////////////////////////
    //                      Question 2                      //
    //////////////////////////////////////////////////////////

    string compile_commandQ2 = "mpiCC -o Tut8_Q2 Tut8_Q2.cpp";
    
    int compile_status_2 = system(compile_commandQ2.c_str());

    if (compile_status_2 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }


    string run_command_2 = "mpiexec -np 1  ./Tut8_Q2";
    cout<<"**********************************************"<<endl;
    cout<<"Running Question 2 with 1 processor"<<endl;
    cout<<"**********************************************"<<endl;
    int run_status_2 = system(run_command_2.c_str());

    if (run_status_2 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }


    //////////////////////////////////////////////////////////
    //                      Question 3                     //
    //////////////////////////////////////////////////////////

    string compile_commandQ3 = "mpiCC -o Tut8_Q3 Tut8_Q3.cpp";
    
    int compile_status_3 = system(compile_commandQ3.c_str());

    if (compile_status_3 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }


    string run_command_3 = "mpiexec -np 6  ./Tut8_Q3";
    cout<<"**********************************************"<<endl;
    cout<<"Running Question 3 with 1 processor"<<endl;
    cout<<"**********************************************"<<endl;
    int run_status_3 = system(run_command_3.c_str());

    if (run_status_3 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }



    //////////////////////////////////////////////////////////
    //                      Question 4                      //
    //////////////////////////////////////////////////////////

    string compile_commandQ4 = "mpiCC -o Tut8_Q4 Tut8_Q4.cpp";
    
    int compile_status_4 = system(compile_commandQ4.c_str());

    if (compile_status_4 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 1;
    }


    string run_command_4 = "mpiexec -np 6  ./Tut8_Q4";
    cout<<"**********************************************"<<endl;
    cout<<"Running Question 4 with 1 processor"<<endl;
    cout<<"**********************************************"<<endl;
    int run_status_4 = system(run_command_4.c_str());

    if (run_status_4 != 0) 
    {
        cerr << "Error running MPI program." << endl;
        return 1;
    }


*/

    //////////////////////////////////////////////////////////
    //                      Question 5                      //
    //////////////////////////////////////////////////////////
    
    string compile_commandQ5 = "mpiCC -o Tut8_Q5 Tut8_Q5.cpp";

    int compile_status_5 = system(compile_commandQ5.c_str());

    if (compile_status_5 != 0) 
    {
        cerr << "Error compiling MPI program." << endl;
        return 5;
    }
    
    vector<string> num_procs_5 = {"2","4"};
    vector<string> ns = {"32","256"};
    for(auto procs : num_procs_5)
    {
        for(auto n : ns)
        {
            string run_command_5 = "mpiexec -np " + procs + " ./Tut8_Q5 "+n;
            cout<<"**********************************************"<<endl;
            cout<<"Running Question 5 with "<<procs<<" number of processors"<<endl;
            cout<<"**********************************************"<<endl;
            int run_status_5 = system(run_command_5.c_str());

            if (run_status_5 != 0) 
            {
                cerr << "Error running MPI program." << endl;
                return 5;
            }
        }
    }
    

    return 0;
}
