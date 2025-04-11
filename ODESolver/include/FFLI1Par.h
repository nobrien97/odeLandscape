#include <vector>
#include <memory>
#include "odePar.h"
#pragma once
class FFLI1Par : public ODEPar
{
private:
    double _AUC = 6.17; // default value when all parameters are 1
public:
    /*FFLI1Par(std::vector<double> traits, std::vector<double> pars);*/
    FFLI1Par();

    std::vector<double> SolveODE() override;

    const static int numPars = 9;
    const static int numTraits = 3;
 
    // Molecular components
    const double& AUC() const { return _AUC; }
    const double& aY() const { return _pars[0]; }
    const double& bY() const { return _pars[1]; }
    const double& KY() const { return _pars[2]; }
    const double& aZ() const { return _pars[3]; }
    const double& bZ() const { return _pars[4]; }
    const double& KXZ() const { return _pars[5]; }
    const double& base() const { return _pars[6]; }
    const double& n() const { return _pars[7]; }
    const double& XMult() const { return _pars[8]; }

    double& AUC() { return _AUC; }
    double& aY() { return _pars[0]; }
    double& bY() { return _pars[1]; }
    double& KY() { return _pars[2]; }
    double& aZ() { return _pars[3]; }
    double& bZ() { return _pars[4]; }
    double& KXZ() { return _pars[5]; }
    double& base() { return _pars[6]; }
    double& n() { return _pars[7]; }
    double& XMult() { return _pars[8]; }

    // Traits
    const double& TimeToHalfMaxExpression() const { return _solutionTraits[0]; }
    const double& MaxExpression() const {return _solutionTraits[1]; }
    const double& TimeAboveHalfMaxExpression() const { return _solutionTraits[2]; }

    const inline void SetTimeToHalfMaxExpression(double value) { _solutionTraits[0] = value; }
    const inline void SetMaxExpression(double value) { _solutionTraits[1] = value; }
    const inline void SetTimeAboveHalfMaxExpression(double value) { _solutionTraits[2] = value; }


};
