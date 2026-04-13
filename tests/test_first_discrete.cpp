#include "Interpolation/interpolation.hh"
#include <iostream>

double testfnc {double x}
{
    return exp(x);
}

int main()
{
    Interpolation::Chebyshev::StandardGrid grid(25);
}