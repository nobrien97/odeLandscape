#include <vector>
#include <memory>
#include "ascent/Ascent.h"
#define STATS_ENABLE_EIGEN_WRAPPERS
#define STATS_ENABLE_STDVEC_WRAPPERS
#include "stats.hpp"
#pragma once
class ODEPar
{
protected:
    double _AUC = 1.0; // default value when all parameters are 1
    std::vector<double> _pars;
    std::vector<double> _solutionTraits;

public:
    ODEPar() {};
    ODEPar(int pars);
    /*ODEPar(int numPar, std::vector<double> pars);*/
    /*ODEPar(int numPar, int numTrait, std::vector<double> traits, std::vector<double> pars); */
    ~ODEPar() = default;

    enum motif_enum {
		NAR,
		PAR,
		FFLC1,
		FFLI1,
		FFBH,
        none
	};

    static motif_enum HashMotifString(std::string motif);

    double static calculateFitness(std::vector<double> pheno, std::vector<double> width, std::vector<double> optimum);
    std::string printPars(std::vector<double> width, std::vector<double> fitnessOptimum, char const *delim);

    void operator ++ (){  
        count++;  
    }  

    bool operator==(const ODEPar rhs) const
    {
        if (numPars != rhs.numPars)
        {
            return false;
        }

        size_t sum = 0;
        for (size_t i = 0; i < numPars; ++i)
        {
            sum += (_pars[i] == rhs._pars[i] ? 1 : 0);
        }
        return sum == numPars;
    }

    static std::unique_ptr<ODEPar> MakeODEPtr(motif_enum motifType);
    static std::unique_ptr<ODEPar> MakeODEPtr(motif_enum motifType, const ODEPar &initialODEPar);

    //static ODEPar* MakeODEPar

    bool Compare(const ODEPar rhs); 

    virtual std::vector<double> SolveODE();

    static double AUC(const double &h, const double &a, const double &b);

    const size_t numPars = 0;
    const size_t numTraits = 0;
    uint count = 1; // Count the number of instances in the population that this exists: needs to be reset to 1 every generation! 
    
    // Set a given value
    void setParValue(int i, double val);

    // Set all values
    void setParValue(std::vector<double> vals, bool firstElementIsAUC = false);

    // Get an ODEPar from a vector of ODEPars
    static double getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec, bool incrementCount = false);

    // Get all the values from an ODEPar
    std::vector<double> getPars(bool returnAUC = true);

    void setAUC(double val) { _AUC = val; }

    double getParValue(int i);

    static std::vector<double> CalcSteadyState(const asc::Recorder &solution, const double& startTime, const double& stopTime, const int &solutionIndex = 0);

    static double CalcDelayTime(const asc::Recorder &solution, const double &startTime, const double &stopTime, const int &solutionIndex, const double& aZ = 1.0, const double& baseline = 0.0);

    static std::vector<double> CalcMaxExpression(const asc::Recorder &solution, const float& startTime, const int &solutionIndex);

    static double CalcTimeAboveThreshold(const asc::Recorder &solution, const double &threshold, const int &solutionIndex);

    static std::vector<double> CalcSecondSteadyState(const asc::Recorder &solution, const double& prevSteadyState, const double& prevSteadyStateTime, const int &solutionIndex);

    // Function for simple baseline regulation when t < Xstart
    static inline double SimpleRegulation(const float& t, const float& aZ, const float &baseline) {
        return (baseline / aZ) * (1 - std::exp(-aZ*t));
    };

    static inline double Interpolate(double x1, double y1, double x2, double y2, double y_target) {
        return x1 + (y_target - y1) * (x2 - x1) / (y2 - y1);
    };
    
    inline std::vector<double> GetTraits() { return _solutionTraits; } 

    inline void SetTraits(std::vector<double> values) 
    {
        if (_solutionTraits.size() != values.size())
        {
            return;
        }

        for (int i = 0; i < _solutionTraits.size(); ++i)
        {
            _solutionTraits[i] = values[i];
        }
    }

};
