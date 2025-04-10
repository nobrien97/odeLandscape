#include "FFBHPar.h"

FFBHPar::FFBHPar(std::vector<double> traits, std::vector<double> pars) : ODEPar(numPars, numTraits, traits, pars)
{
    _pars.resize(numPars, 1.0);
    setParValue(pars, false);
    _solutionTraits.resize(numTraits);
    SetTraits(traits);
    //_pars[8] = 0.0; // set baseline to 0 to start
}

FFBHPar::FFBHPar() : ODEPar(numPars) 
{
    _pars.resize(numPars, 1.0);
    _pars[8] = 0.0; // set baseline to 0 to start
    _solutionTraits.resize(numTraits);
}


std::vector<double> FFBHPar::SolveODE()
{
    // Initialise output
    //std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;

    // Starting X
	int X = 0;

    // Precompute K^n parameters (reduce pow calls)
    const double KYn = pow(KY(), n());
    const double KZXn = pow(KZX(), n());
    const double KXZn = pow(KXZ(), n());

    // FFBH system
	auto FFBHDerivative = [this, &Xstart, &Xstop, &X, &KYn, &KZXn, &KXZn](const asc::state_t &curState, asc::state_t &nextState, double t)
	{
        // # Manually set X
        //     X <- XMult * (t > Xstart && t <= Xstop)
        //     X <- X + XH

        //     # Hill function component of X, XH
        //     dXH <- ( Z^Hilln / (KZX^Hilln + Z^Hilln) ) - aX*XH
            
        //     # Update X
        //     X <- X + dXH
            
        //     dY <- base * X + bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
        //     dZ <- base * X + bZ *  ((X * Y)^Hilln)/( (KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln) ) - aZ*Z
        
        double Xnew = X * XMult();
        Xnew += curState[0];

        double baseline = std::max(base() - 1.0, 0.0); // Adjust baseline so it is relative to the default value, 0 (instead of 1)

        // Hill function/feedback component of X
        nextState[0] = ( pow(curState[2], n()) / (KZXn + pow(curState[2], n())) ) - aX() * curState[0];
        
        // Update X with the Hill function X component
        Xnew += nextState[0];

        // Y
		nextState[1] = bY() * pow(Xnew, n()) / ( KYn + pow(Xnew, n()) ) - aY() * curState[1];
		
        // Z
        nextState[2] = baseline + bZ() * pow(Xnew * curState[1], n()) / ( ( KXZn + pow(Xnew, n()) ) * ( KYn + pow(curState[1], n()) ) ) - aZ() * curState[2];    
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

    // Calculate traits
    std::vector<double> maxExp = ODEPar::CalcMaxExpression(recorder, Xstart, 4);
    std::vector<double> secondSteadyState = ODEPar::CalcSecondSteadyState(recorder, maxExp[0], Xstop, 4);
    

    SetTimeToHalfMaxExpression(maxExp[1]);
    SetMaxExpression(maxExp[0]);
    SetSecondSteadyState(secondSteadyState[1]);
    SetTimeToSecondSteadyState(secondSteadyState[0]);

	// Calculate AUC
	// double z = 0;
	// #pragma omp simd reduction(+:z)
	// for (uint i = 0; i < recorder.history.size()-1; ++i)
	// {
	// 	z += ODEPar::AUC(0.1, (double)recorder.history[i][4], (double)recorder.history[i + 1][4]);
	// }
	
	// Check that z is > 0, set AUC, return
	// z = (z >= 0) ? z : 0.0;
    // result[0] = z;
    // setAUC(z);
    return {TimeToMaxExpression(), MaxExpression(), TimeToSecondSteadyState(), SecondSteadyState()};
}
