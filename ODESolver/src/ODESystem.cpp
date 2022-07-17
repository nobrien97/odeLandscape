#include "ODESystem.h"
#include <numeric>



double ODESystem::AUC(const double &h, const double &a, const double &b)
{
    return ((a + b) * 0.5) * h;
};


double ODESystem::calculatePhenotype()
{
    auto NAR = [this](const asc::state_t &val, asc::state_t &dxdt, double t)
    {
      // dX <- if (t == Xstart) 1; if (t == Xstop) -1; else 0;
      // dZ <- bZ * X^nXZ/(KXZ^nXZ + X^nXZ) * (KZ^nZ/(Kz^nZ+Z^nZ)) - aZ*Z
      dxdt[0] = _pars.bZ() * pow(X, nXZ) / (pow(_pars.KXZ(), nXZ) + pow(X, nXZ)) * (pow(_pars.KZ(), nZ)/(pow(_pars.KZ(), nZ) + pow(val[0], nZ))) - _pars.aZ() * val[0];
    };

    asc::state_t state = { 0.0 };
    double t = 0.0;
    double dt = 0.1;
    double t_end = 10.0;

    asc::RK4 integrator;
    asc::Recorder recorder;

    // Solve
    while (t < t_end)
    {
        // Add a small epsilon to X to get around floating point imprecision
        X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
        recorder({ t, (asc::value_t)X, state[0] });
        integrator(NAR, state, t, dt);
    }

    // Calculate AUC
    std::vector<double> z_auc = std::vector<double>(recorder.history.size());
    for (uint i = 0; i < recorder.history.size()-1; ++i)
    {
        z_auc[i] = AUC(0.1, (double)recorder.history[i][2], (double)recorder.history[i + 1][2]);
    }

    double z = std::accumulate(z_auc.begin(), z_auc.end(), 0.0);
    // Make sure we're above 0.0
    z = (z >= 0) ? z : 0.0;
    this->_pars.setAUC(z);
    return z;
}

std::string ODESystem::printPars(char* const delim)
{
    std::vector<double> parameters = _pars.getPars(); 
    std::string result;
    for (int i = 0; i < size(parameters); ++i)
    {
        if ((size(parameters) - 1) == i) 
        {
            result.append(std::to_string(parameters[i]));
            break;
        }
        result.append(std::to_string(parameters[i]) + delim);
    }

    return result;
}


ODESystem::ODESystem(/* args */)
{

}

ODESystem::~ODESystem()
{
}
