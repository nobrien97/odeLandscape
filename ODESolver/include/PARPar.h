#include <vector>
#include <memory>
#include "odePar.h"
#pragma once
class PARPar : public ODEPar
{
private:
    double _AUC = 6.17; // default value when all parameters are 1
public:
    /*PARPar(std::vector<double> traits, std::vector<double> pars);*/
    PARPar();

    std::vector<double> SolveODE() override;

    const static int numPars = 7;
    const static int numTraits = 2;
 
    // Molecular components
    const double& AUC() const { return _AUC; }
    const double& aZ() const { return _pars[0]; }
    const double& bZ() const { return _pars[1]; }
    const double& KZ() const { return _pars[2]; }
    const double& KXZ() const { return _pars[3]; }
    const double& base() const { return _pars[4]; }
    const double& n() const { return _pars[5]; }
    const double& XMult() const { return _pars[6]; }

    double& AUC() { return _AUC; }
    double& aZ() { return _pars[0]; }
    double& bZ() { return _pars[1]; }
    double& KZ() { return _pars[2]; }
    double& KXZ() { return _pars[3]; }
    double& base() { return _pars[4]; }
    double& n() { return _pars[5]; }
    double& XMult() { return _pars[6]; }

    // Traits
    const double& ResponseTime() const { return _solutionTraits[0]; }
    //const double& ResponseDelay() const {return _solutionTraits[1]; } // NOTE: ResponseDelay requires very specific conditions to arise, secondary function, removing for now
    const double& SteadyState() const { return _solutionTraits[1]; }

    const inline void SetResponseTime(double value) { _solutionTraits[0] = value; }
    //const inline void SetResponseDelay(double value) { _solutionTraits[1] = value; }
    const inline void SetSteadyState(double value) { _solutionTraits[1] = value; }

};
