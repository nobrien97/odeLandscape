#include "FFBHPar.h"

FFBHPar::FFBHPar(double AUC, std::vector<double> pars) : ODEPar(numPars, pars)
{
    _AUC = AUC;
    _pars.resize(numPars, 1.0);
    _pars[9] = 0.0; // set baseline to 0 to start
}

FFBHPar::FFBHPar() : ODEPar(numPars) 
{
    _pars.resize(numPars, 1.0);
    _pars[9] = 0.0; // set baseline to 0 to start

}


std::vector<double> FFBHPar::SolveODE()
{
    // Initialise output
    std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;

    // Starting X
	int X = 0;

    // FFBH system
	auto FFBHDerivative = [this, &Xstart, &Xstop, &X](const asc::state_t &curState, asc::state_t &nextState, double t)
	{

        // # Manually set X
        // X <- XMult * (t > Xstart && t <= Xstop)
        // X <- X + XH
        
        // # Hill function component of X, XH
        // dXH <- ( Z^Hilln / (KZX^Hilln + Z^Hilln) ) - aX*XH
        
        // # Update X
        // X <- X + dXH
        
        // dY <- base * X + bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
        // dZ <- base * X + bZ *  ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
        // Update X with the Hill function X component
        double Xnew = XMult() * X;
        Xnew += curState[0];

        // Hill function/feedback component of X
        nextState[0] = ( pow(curState[2], n()) / (pow(KZX(), n()) + pow(curState[2], n())) ) - aX() * curState[0];

        // Update Xnew with dXH
        Xnew += nextState[0];
        
        // Y
		nextState[1] = base() * Xnew + bY() * pow(Xnew, n()) / ( pow(KY(), n()) + pow(Xnew, n()) ) - aY() * curState[1];
		
        // Z
        nextState[2] = base() * Xnew + bZ() * pow(Xnew * curState[1], n()) / ( ( pow(KXZ(), n()) + pow(Xnew, n()) ) * ( pow(KY(), n()) + pow(curState[1], n()) ) ) - aZ() * curState[2];    
    };

	// Set up the initial state
	asc::state_t state = { 0.0, 0.0, 0.0 };
	double t = 0.0;
	double dt = 0.1;
	double t_end = 10.0;
	asc::RK4 integrator;
	asc::Recorder recorder;

	while (t < t_end)
	{
		// Add a small epsilon to get around t floating point inaccuracy
		X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
		recorder({t, (asc::value_t)X, state[0], state[1], state[2]});
		integrator(FFBHDerivative, state, t, dt);
	}
	// Calculate AUC
	double z = 0;
	#pragma omp simd reduction(+:z)
	for (uint i = 0; i < recorder.history.size()-1; ++i)
	{
		z += ODEPar::AUC(0.1, (double)recorder.history[i][4], (double)recorder.history[i + 1][4]);
	}
	
	// Check that z is > 0, set AUC, return
	z = (z >= 0) ? z : 0.0;
    result[0] = z;
    setAUC(z);
    return result;
}
