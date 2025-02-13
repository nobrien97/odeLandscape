#include "FFLC1Par.h"

FFLC1Par::FFLC1Par(double AUC, std::vector<double> pars) : ODEPar(numPars, pars)
{
    _AUC = AUC;
    _pars.resize(numPars, 1.0);
    _pars[6] = 0.01; // constitutive promoter
}

FFLC1Par::FFLC1Par() : ODEPar(numPars) 
{
	_pars.resize(numPars, 1.0);
    _pars[6] = 0.01; // constitutive promoter
}

std::vector<double> FFLC1Par::SolveODE()
{
    // Initialise output
    std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;

    // Starting X and Y
	int X = 0;

    // FFL C1 system
	auto FFLC1Derivative = [this, &Xstart, &Xstop, &X](const asc::state_t &curState, asc::state_t &nextState, double t)
	{
        // X <- XMult * (t > Xstart && t <= Xstop)
        // dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
        // dZ <- bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
        
        double Xnew = X * XMult();

		nextState[0] = base() * Xnew + bY() * pow(Xnew, n()) / (pow(KY(), n()) + pow(Xnew, n())) - aY() * curState[0];
		nextState[1] = base() * Xnew + bZ() * pow(Xnew * curState[0], n()) / ( (pow(KXZ(), n()) + pow(Xnew, n())) * (pow(KY(), n()) + pow(curState[0], n())) ) - aZ() * curState[1];    
    };

	// Set up the initial state
	asc::state_t state = { 0.0, 0.0 };
	double t = 0.0;
	double dt = 0.1;
	double t_end = 10.0;
	asc::RK4 integrator;
	asc::Recorder recorder;

	while (t < t_end)
	{
		// Add a small epsilon to get around t floating point inaccuracy
		X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
		recorder({t, (asc::value_t)X, state[0], state[1]});
		integrator(FFLC1Derivative, state, t, dt);
	}
	// Calculate AUC
	double z = 0;
	#pragma omp simd reduction(+:z)
	for (uint i = 0; i < recorder.history.size()-1; ++i)
	{
		z += ODEPar::AUC(0.1, (double)recorder.history[i][3], (double)recorder.history[i + 1][3]);
	}
	
	// Check that z is > 0, set AUC, return
	z = (z >= 0) ? z : 0.0;
    result[0] = z;
    setAUC(z);
    return result;
}