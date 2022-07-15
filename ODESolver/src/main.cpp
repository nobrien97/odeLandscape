#include "ascent/Ascent.h"
#include "ODESystem.h"
#include "odePar.h"
#include <vector>
#include <iostream>

int main(int argc, char* argv[])
{
    //TODO: Fix this, get command line input working - make it work on a file and extract the values from there, operate over a range (e.g. specify aZ min to aZ max)
    ODEPar parameters(-1.0, argv[1], argv[2], argv[3], argv[4]);

    ODESystem ODE(parameters);

    // Calculate phenotype, write to file
    ODE.calculatePhenotype();

    return 0;
}