#include <iostream>
#include "gnuplot-iostream.h"

int main() 
{
    Gnuplot gp;
    gp << "plot sin(x)*tan(x)\n";
    return 0;
}
