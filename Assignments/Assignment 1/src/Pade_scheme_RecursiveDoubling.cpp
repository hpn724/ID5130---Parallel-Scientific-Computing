#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>

using namespace std;


//Function to show matrix
void show_matrix(double** A, int m, int n) 
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}


//Function to show vector
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

//Function to get the actual values of the defined function
double* get_func_values(double a, double b, int n, double (*func)(double))
{
    double h = (b-a)/(double(n-1));
    double* f = new double[n];
    int i;
    for(i=0;i<n;i++)
    {
        f[i]=func(a+i*h);
    }
    return f;

}

//defined function to be used
double function(double x) 
{
    double f=sin(5*x);
    return f;
}

//actual derivative of the funtion
double actual_der(double x)
{
    return 5*cos(5*x);
}

//actual derivative values of the function
double* get_actual_sol(double a, double b, int n, double (*f)(double))
{
    double* y;
    double h = (b-a)/(double(n-1));
    y = new double[n];
    int i;
    for(i=0;i<n;i++)
    {
        y[i]=f(a+i*h);
    }
    return y;
}

//function to initialize the pade matrix and rhs vector
void initialize_pade_scheme_variables(double**& A, double*& rhs, double* f, int n, double a, double b) 
{
    double h = (b - a)/(double(n-1));

    rhs = new double[n ];

    int i,j;
    A = new double*[n + 1];

    for (i = 0; i <n; i++) {
        A[i] = new double[n + 1];
        for (j = 0; j <n; j++) {
            A[i][j] = 0;
        }
    }

    for (i = 0; i < n; i++) 
    {
        if (i == 0) 
        {
            rhs[i] = (1 / h) * (-2.5 * f[0] + 2 * f[1] + 0.5 * f[2]);
        } 
        else if (i == n-1) 
        {
            rhs[i] = (1 / h) * (2.5 * f[n-1] - 2 * f[n-2] - 0.5 * f[n-3]);
        } 
        else 
        {
            rhs[i] = (3 / h) * (f[i+1] - f[i-1]);
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
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



//Recursive doubling algorithm which takes the number of threads as a parameter
void recursive_doubling(int n,double** A, double* yik, double* &result, int thread,double* time_taken)
{
    double timer=0;
    double* aik = new double[n];
    double* bik = new double[n];
    double* cik = new double[n];
    double* aik1 = new double[n];
    double* bik1 = new double[n];
    double* cik1 = new double[n];
    double* yik1 = new double[n];
    double* alphaki = new double[n];
    double* betaki = new double[n];

    int l1;
    int l2;
    int i;

    timer=omp_get_wtime();
    
    #pragma omp parallel for num_threads(thread) private(i) shared(aik,bik,cik,n,A)
    for(i=0;i<n-1;i++)
    {
        aik[i+1] = A[i+1][i];
        bik[i] = A[i][i];
        cik[i] = A[i][i+1];
    }

    bik[n-1]=A[n-1][n-1];

    
    //Start of parallelisation
    #pragma omp parallel num_threads(thread) shared(A,aik,bik,cik,yik,aik1,bik1,cik1,alphaki, betaki,n) private(i, l1, l2)
    {
        for(int k=0;k<=ceil(log2(n));k++)
        {
            l1 = pow(2,k-1);
            l2 = pow(2,k);

            #pragma omp for
                for(i=0;i<n;i++)
                {
                    bik1[i]=bik[i];
                    yik1[i]=yik[i];

                    if(i <= n - l1 - 1) 
                    {
                        betaki[i] = -cik[i]/bik[i+l1];

                        bik1[i] += betaki[i]*aik[i+l1];
                        yik1[i] += betaki[i]*yik[i+l1];
                    }
                    else 
                    {
                        betaki[i] = 0;
                    }

                    if (i >= l1) 
                    {
                        alphaki[i] = -aik[i]/bik[i-l1];

                        bik1[i] += alphaki[i]*cik[i-l1];
                        yik1[i] += alphaki[i]*yik[i-l1];
                    }
                    else 
                    {
                        alphaki[i] = 0;
                    }

                    if (i <= n - l2 - 1) 
                    {
                        cik1[i] = betaki[i]*cik[i+l1];
                    }
                    else 
                    {
                        cik1[i] = 0;
                    }

                    if (i >= l2) 
                    {
                        aik1[i] = alphaki[i]*aik[i-l1];
                    }
                    else 
                    {
                        aik1[i] = 0;
                    }
                }
            
            if (k <= ceil(log2(n)) - 1) 
            {
                #pragma omp for
                    for (i = 0; i < n; i++) 
                    {
                        aik[i] = aik1[i];
                        bik[i] = bik1[i];
                        cik[i] = cik1[i];
                        yik[i] = yik1[i];
                    }
            
            }
        }
    }
    //end of parallel block
    

    result = new double[n];

    #pragma omp parallel for num_threads(thread) private(i) shared(result,yik1,bik1,n)
    for(int i=0;i<n;i++)
    {
        result[i]=yik1[i]/bik1[i];
    }
    
    timer = omp_get_wtime() - timer;

    *time_taken=timer;
}


//function to get the RD data as a csv file
void get_plot_data(double a, double b, int n, double* y1, double* y2)
{
    double* x;
    double h = (b-a)/(double(n-1));
    x = new double[n];
    for(int i=0;i<n;i++)
    {
        x[i] = a+i*h;
    }

    std::ofstream dataFile("Pade_scheme_RecursiveDoubling_data.csv");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        if(i==0)
        {
            dataFile << "x,Actual value,Recursive Doubling" << endl;
        }
        dataFile << x[i] << "," << y1[i] << "," << y2[i] << endl;
    }
    dataFile.close();
}

//function to get the time vs threads data as a csv file
void get_time_plot_data(double* thread, double* time,  int n)
{
    std::ofstream dataFile("PadeScheme_RecursiveDoubling_Thread_Time.csv");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        if(i==0)
        {
            dataFile << "Thread Count,Time" << endl;
        }
        dataFile << thread[i] << "," << time[i] << endl;
    }
    dataFile.close();
}

int main()
{
    double a = 0.0;
    double b = 3.0;
    int n = 100;
    double** A;
    double* rhs;

    
    
    int thread = 2;

    double* yik = get_func_values(a,b,n,&function);
    double* dy = get_actual_sol(a, b, n, &actual_der);

    initialize_pade_scheme_variables(A, rhs, yik, n, a, b);

    double* result;
    double time_taken;
    recursive_doubling(n, A, rhs, result, thread,&time_taken);
    
    cout<<"Time taken for "<<thread<<" threads : "<<time_taken<<endl;

    get_plot_data(a, b, n, dy, result);


    //To get time dependency on thread count

    double threads[]={2,4,8};
    double* time_array;
    int N=10000;
    time_array = new double[3];
    
    double** AA;
    double* RHS;


    double* Y = get_func_values(a,b,N,&function);

    time_array = new double[3];

    initialize_pade_scheme_variables(AA, RHS, yik, N, a, b);

    cout<<"For n = 1000 :"<<endl;
    for(int i=0;i<3;i++)
    {
        double* temp_result;
        double temp_time=0;
        recursive_doubling(n, AA, RHS, temp_result, threads[i],&temp_time);
        cout<<"Time taken for "<<threads[i]<<" threads : "<<temp_time<<endl;
        time_array[i]=temp_time;
        delete [] temp_result;

    }
    
    get_time_plot_data(threads, time_array, 3);
        
        delete [] RHS;
        for(int i=0;i<N;i++)
        {
            delete [] AA[i];
        }
        delete [] AA;
    delete [] time_array;
}