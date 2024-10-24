#include <vector>
#include <memory>
#include "ascent/Ascent.h"
#pragma once
class ODEPar
{
protected:
    double _AUC = 1.0; // default value when all parameters are 1
    std::vector<double> _pars;
public:
    ODEPar() {};
    ODEPar(int pars);
    ODEPar(int numPar, std::vector<double> pars);
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

    double static calculateFitness(double pheno, double width, double optimum);
    std::string printPars(double width, double fitnessOptimum, char const *delim);

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
    unsigned int count = 1; // Count the number of instances in the population that this exists: needs to be reset to 1 every generation! 
    
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

};
