#include<iostream>
#include<math.h>
#include <vector>

#include<omp.h>

using namespace std;

double function(double x)
{
    return sin(x)/(2*pow(x,3));
}

double trapezoidal_integral(double a, double b, double(*func)(double), int n)
{
    double dx=(b-a)/n;
    double s1=func(a);
    
    double sn=func(b);
    double sab=0;
    for(double i=a+dx;i<b;i+=dx)
    {
        sab+=2*func(i);
        
    }
    return dx*0.5*(s1+sab+sn);
}

void trap_parallel_critical(double a, double b, double(*func)(double), int n, double *value)
{
    double dx=(b-a)/n;
    int thread_count=omp_get_num_threads();
    int thread_rank=omp_get_thread_num();
    int local_n=n/thread_count;
    double local_a=a+thread_rank*local_n*dx;
    double local_b=local_a+local_n*dx;

    double local_ints=(func(a)+func(b))/2;
    for(double i=local_a+dx;i<local_b;i+=dx)
    {
        local_ints+=func(i);
        
    }
    local_ints=local_ints*dx;
    
    #pragma omp critical
    {
        *value+=local_ints;
    }

    return;
    
}

double trap_parallel_for_red(double a, double b, double(*func)(double), int n,int number_thread)
{
    double dx=(b-a)/n;
    double s1=func(a);
    
    double sn=func(b);
    double sab=0;
    #pragma omp parallel for num_threads(number_thread) reduction(+ :sab)
    for(int i=1;i<n;i++)
    {
        double x = a + i*dx;
        sab+=2*func(x);

    }

    return dx*0.5*(s1+sab+sn);
}


int main()
{
    double a,b;
    a=1;
    b=M_PI;
    int n=32;

    cout<<endl;

    cout<<"############# Normal for loop to find the area #############"<<endl;

    cout<<endl;
    

    double trap_int=trapezoidal_integral(a,b,&function,n);


    cout<<"The trapezoidal integral value is : "<<trap_int<<endl;
    
    vector<int> num_thread={2,4,8};
    
    cout<<endl;

    cout<<"############# Using omp parallel critical to find the area #############"<<endl;

    cout<<endl;

    for(int i=0;i<num_thread.size();i++)
    {
        double trap_int_crit=0.0;
        #pragma omp parallel num_threads(num_thread[i])
        trap_parallel_critical(a, b, &function, n,&trap_int_crit);
        cout<<"The trapezoidal integral parallel critical for number of threads : " <<num_thread[i] <<" = "<<trap_int_crit<<endl;
    }

    cout<<endl;

    cout<<"############# Using omp parallel for and reduction to find the area #############"<<endl;

    cout<<endl;

    for(int i=0;i<num_thread.size();i++)
    {
        double trap_int_red_for=0.0;
        trap_int_red_for = trap_parallel_for_red(a, b, &function, n,num_thread[i]);
        cout<<"The trapezoidal integral parallel critical for number of threads : " <<num_thread[i] <<" = "<<trap_int_red_for<<endl;
    }

    


    return 0;
}