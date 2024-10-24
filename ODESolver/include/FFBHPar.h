#include <vector>
#include <memory>
#include "odePar.h"
#pragma once
class FFBHPar : public ODEPar
{
private:
    double _AUC = 2.5353073611315153; // default value when all parameters are 1
public:
    FFBHPar(double AUC, std::vector<double> pars);
    FFBHPar();

    std::vector<double> SolveODE() override;

    const static int numPars = 11;

    const double& AUC() const { return _AUC; }
    const double& aX() const { return _pars[0]; }
    const double& KZX() const { return _pars[1]; }
    const double& aY() const { return _pars[2]; }
    const double& bY() const { return _pars[3]; }
    const double& KY() const { return _pars[4]; }
    const double& aZ() const { return _pars[5]; }
    const double& bZ() const { return _pars[6]; }
    const double& KXZ() const { return _pars[7]; }
    const double& base() const { return _pars[8]; }
    const double& n() const { return _pars[9]; }
    const double& XMult() const { return _pars[10]; }

    double& AUC() { return _AUC; }
    double& aX() { return _pars[0]; }
    double& KZX() { return _pars[1]; }
    double& aY() { return _pars[2]; }
    double& bY() { return _pars[3]; }
    double& KY() { return _pars[4]; }
    double& aZ() { return _pars[5]; }
    double& bZ() { return _pars[6]; }
    double& KXZ() { return _pars[7]; }
    double& base() { return _pars[8]; }
    double& n() { return _pars[9]; }
    double& XMult() { return _pars[10]; }

};
