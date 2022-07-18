#include "odePar.h"
#include "ascent/Ascent.h"

class ODESystem
{
private:
    /* data */
    ODEPar _pars;
    const double Xstart = 1.0;
    const double Xstop = 6.0;
    const double nXZ = 8.0;
    const double nZ = 8.0;
    int X = 0;

public:
    ODESystem(/* args */);
    ODESystem(ODEPar pars) : _pars(pars) {};
    ~ODESystem();

    static double AUC(const double &h, const double &a, const double &b);

    double calculatePhenotype();

    double calculateFitness(double width, double optimum);

    std::string printPars(double width, double fitnessOptimum, char* const delim = ",");
};
