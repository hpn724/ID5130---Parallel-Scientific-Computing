#include<iostream>
#include<math.h>
#include<sstream>
#include<fstream>
#include<string>
#include<gnuplot-iostream.h>


using namespace std;

void show_matrix(double** A, int m, int n) 
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void show_vector(double* A, int n) 
{
    cout << "[ ";
    for (int i = 0; i < n; i++) {
        cout << A[i];
        if (i != n - 1) {
            cout << ", ";
        }
    }
    cout << " ]" << endl;
}


double* get_func_values(double a, double b, int n, double (*func)(double))
{
    double h = (b-a)/(double(n-1));
    double* f = new double[n];
    for(int i=0;i<n;i++)
    {
        
        f[i]=func(a+i*h);
    }
    return f;

}


void initialize_pade_scheme_variables(double**& A, double*& rhs, double* f, int n, double a, double b) 
{
    double h = (b - a)/(double(n-1));

    rhs = new double[n];

    A = new double*[n + 1];
    for (int i = 0; i <n; i++) {
        A[i] = new double[n + 1];
        for (int j = 0; j <n; j++) {
            A[i][j] = 0;
        }
    }

    rhs[0] = (1 / h) * (-2.5 * f[0] + 2 * f[1] + 0.5 * f[2]);
    rhs[n-1] = (1 / h) * (2.5 * f[n-1] - 2 * f[n-2] - 0.5 * f[n-3]);
    for (int i = 1; i < n-1; i++) 
    {
        rhs[i] = (3 / h) * (f[i+1] - f[i-1]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                if (j == 0) {
                    A[i][j] = 1;
                } else if (j == 1) {
                    A[i][j] = 2;
                }
            } else if (i == n-1) {
                if (j == n-1) {
                    A[i][j] = 1;
                } else if (j == n - 2) {
                    A[i][j] = 2;
                }
            } else {
                if (j == i - 1 || j == i + 1) {
                    A[i][j] = 1;
                } else if (i == j) {
                    A[i][j] = 4;
                }
            }
        }
    }
}

void LU_decomposition(double** A, double**& L, double**& U, int n) 
{
    L = new double*[n];
    U = new double*[n];

    for (int i = 0; i < n; i++) {
        L[i] = new double[n];
        U[i] = new double[n];
        for (int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<i)
            {
                L[j][i]=0;
            }
            else
            {
                L[j][i]=A[j][i];
                for(int k=0;k<i;k++)
                {
                    L[j][i]=L[j][i] - L[j][k] * U[k][i];
                }
            }
        }

        for(int j=0;j<n;j++)
        {
            if(j<i)
            {
                U[i][j]=0;
            }
            else if(j==i)
            {
                U[i][j]=1;
            }
            else
            {
                U[i][j]=A[i][j]/L[i][i];
                for(int k=0;k<i;k++)
                {
                    U[i][j]= U[i][j] - ((L[i][k]*U[k][j])/L[i][i]);
                }
            }
        }
    }
}

double* forward_substitution(double** L, double* b, int n)
{
    double temp;
    double* x;
    x=new double[n];
    for(int i=0;i<n;i++)
    {
        x[i]=0;
    }
    x[0]=b[0];
    for(int i=1;i<n;i++)
    {
        double sum=0;
        for(int j=0;j<i;j++)
        {
            sum+=L[i][j]*x[j];
        }
        x[i]=(b[i]-sum)/L[i][i];
    }
    return x;
}

double* backward_substitution(double** U, double* b, int n)
{
    double temp;
    double* x;
    x=new double[n];
    for(int i=0;i<n;i++)
    {
        x[i]=0;
    }
    x[n-1] = b[n-1]/U[n-1][n-1];
    for(int i=n-2;i>=0;i--)
    {
        temp=b[i];
        for(int j=i+1;j<n;j++)
        {
            temp -= U[i][j]*x[j];
        }
        x[i] = temp/U[i][i];
    }
    return x;
}

double func(double x) 
{
    return sin(5*x);
}

double actual_der(double x)
{
    return 5*cos(5*x);
}

double* LU_solve(double** L, double** U, double* b, int n)
{
    double* y = forward_substitution(L, b, n);
    double* x = backward_substitution(U, y, n);
    return x;
}

double* get_actual_sol(double a, double b, int n, double (*f)(double))
{
    double* y;
    double h = (b-a)/(double(n-1));
    y = new double[n];
    for(int i=0;i<n;i++)
    {
        y[i]=f(a+i*h);
    }
    return y;
}

string get_plot_data(double a, double b, int n, double* y1, double* y2)
{
    double* x;
    double h = (b-a)/(double(n-1));
    x = new double[n];
    for(int i=0;i<n;i++)
    {
        x[i] = a+i*h;
    }

    string plot_name = "Pade_scheme_LU_Decomp_plot_data.csv";

    ofstream dataFile(plot_name);
    if (!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return "";
    }
    for (int i = 0; i < n; i++) {
        if(i==0)
        {
            dataFile << "x,Actual value,LU Decomposition values" << std::endl;
        }
        dataFile << x[i] << ","<< y1[i] << ","<< y2[i] << std::endl;
    }
    dataFile.close();
    return plot_name;
}

/*
void plot(const string& filename)
{
    ifstream file(filename);
    if(!file.is_open())
    {
        std::cerr<<"Error: unable to open file : "<<filename<<endl;
        return;
    }

    string header;
    getline(file, header);
    istringstream headerStream(header);
    string columnHeader;
    headerStream >> columnHeader;

    ofstream gnuplotScript("PadeLUdecomp_plot_script.plt");
    gnuplotScript << "set title 'Pade Scheme using LU decomposition'\n";
    gnuplotScript << "set xlabel 'x'\n";
    gnuplotScript << "set ylabel '"<<columnHeader<<"'\n";
    gnuplotScript << "plot from"

    int datasetCount = 0;
    while (headerStream >> columnHeader) 
    {
        if (datasetCount > 0) 
        {
            gnuplotScript << ", ";
        }
        gnuplotScript << "'" << filename << "' using 1:" << datasetCount + 2 << " with linespoints title '" << columnHeader << "'";
        datasetCount++;
    }
    file.close();
    gnuplotScript.close();

    system("gnuplot -persistent PadeLUdecomp_plot_script.plt");

    
}

*/

int main() 
{
    double a = 0.0;
    double b = 3.0;
    int n = 25;
    double** A;
    double* rhs;

    double* y = get_func_values(a,b,n, &func);
    initialize_pade_scheme_variables(A, rhs, y, n, a, b);

    double** L; 
    double** U;

    LU_decomposition(A, L, U, n);

    double* pade_solution = LU_solve(L, U, rhs, n);

    double* dy = get_actual_sol(a, b, n, &actual_der);

    string data = get_plot_data(a, b, n, dy, pade_solution);

    //plot(data);
    
    // Free dynamically allocated memory
    for (int i = 0; i < n; ++i)
    {
        delete[] A[i];
        delete[] U[i];
        delete[] L[i];
    }
    delete[] A;
    delete[] U;
    delete[] L;
    delete[] rhs;
    delete[] pade_solution;
    delete[] y;

    return 0;


}