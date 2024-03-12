#include <iostream>
#include<stdio.h>


int main() {
    std::cout << "Hello Easy C++ project!" << std::endl;
    #pragma omp parallel for
    {
        for(int i=0;i<10;i++)
        {
            printf(" line %d \n",i);
        }
    }
    return 0;
}