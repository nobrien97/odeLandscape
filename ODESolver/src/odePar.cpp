#include "odePar.h"

ODEPar::ODEPar(/* args */)
{
}

ODEPar::~ODEPar()
{
}

void ODEPar::setParValue(size_t i, double val)
{
    switch (i)
    {
    case 0:
        _AUC = val;
        break;
    case 1:
        _aZ = val;
        break;
    case 2:
        _bZ = val;
        break;
    case 3:
        _KZ = val;
        break;
    case 4:
        _KXZ = val;
    default:
        break;
    }

    return;
}

void ODEPar::setParValue(std::vector<double> vals)
{
    if (vals.size() != 5)
        return;
    
    for (int i = 0; i < 5; ++i)
    {
        setParValue(i, vals[i]);
    }
}

std::vector<double> ODEPar::getPars()
{
    return std::vector<double>({ _AUC, _aZ, _bZ, _KZ, _KXZ });
} 

// Get an ODEPar from a vector of ODEPars
double ODEPar::getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec, bool incrementCount)
{
    for (auto& ODE : vec)
    {
        if (target == *ODE.get())
        {
            if (incrementCount) 
                ++*ODE;

            return ODE->AUC();
        }
    }
    return 0;
} 
