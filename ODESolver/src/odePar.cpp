#include "odePar.h"
#include "NARPar.h"
#include "PARPar.h"
#include "FFLC1Par.h"
#include "FFLI1Par.h"
#include "FFBHPar.h"

#define MAX_EXP 10000
#define MAX_TIME 10.0

ODEPar::motif_enum ODEPar::HashMotifString(std::string motif)
{
    if (motif == "NAR") return motif_enum::NAR;
    if (motif == "PAR") return motif_enum::PAR;
    if (motif == "FFLC1") return motif_enum::FFLC1;
    if (motif == "FFLI1") return motif_enum::FFLI1;
    if (motif == "FFBH") return motif_enum::FFBH;
    
    return motif_enum::none;
}

std::string ODEPar::printPars(std::vector<double> width, std::vector<double> fitnessOptimum, char const *delim = ",")
{
    // Calculate fitness to print
    std::vector<double> traits = GetTraits();

    double fitness = calculateFitness(traits, width, fitnessOptimum);
    std::vector<double> parameters = getPars(false);

    // Insert trait values and fitness to the parameters array
    parameters.insert(parameters.begin(), traits.begin(), traits.end());
    parameters.insert(parameters.begin(), fitness); 
    std::string result;
    for (int i = 0; i < parameters.size(); ++i)
    {
        if ((parameters.size() - 1) == i) 
        {
            result.append(std::to_string(parameters[i]));
            break;
        }
        result.append(std::to_string(parameters[i]) + delim);
    }
    return result;
}

double ODEPar::calculateFitness(std::vector<double> pheno, std::vector<double> width, std::vector<double> optimum)
{
    // Convert width to v/cov matrix
    // TODO: This might error because DiagonalMatrix has EigenBase as a base class, not MatrixBase
    // Seems like some functions on MatrixBase are used in stats::dmvnorm?
    Eigen::VectorXd width_xd = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(width.data(), width.size());
    Eigen::VectorXd pheno_xd = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(pheno.data(), pheno.size());
    Eigen::VectorXd optimum_xd = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(optimum.data(), optimum.size());

    Eigen::MatrixXd Sigma = width_xd.asDiagonal().toDenseMatrix();

    double fitness = stats::dmvnorm(pheno_xd, optimum_xd, Sigma);
    double normConst = stats::dmvnorm(optimum_xd, optimum_xd, Sigma);
    return fitness / normConst;
}


std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType)
{
        {
        switch (motifType)
        {
            case NAR:
                return std::make_unique<NARPar>();
                break;
            case PAR:
                return std::make_unique<PARPar>();
                break;
            case FFLC1:
                return std::make_unique<FFLC1Par>();
                break;
            case FFLI1:
                return std::make_unique<FFLI1Par>();
                break;
            case FFBH:
                return std::make_unique<FFBHPar>();
                break;
            default:
                return nullptr;
                break;
        }
    }
}
/*
std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType, const ODEPar &initialODEPar)
{
    switch (motifType)
    {
        case NAR:
            return std::make_unique<NARPar>(7, 2, initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case PAR:
            return std::make_unique<PARPar>(7, 2, initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFLC1:
            return std::make_unique<FFLC1Par>(9, 3, initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFLI1:
            return std::make_unique<FFLI1Par>(9, 3, initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFBH:
            return std::make_unique<FFBHPar>(11, 4, initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        default:
            return nullptr;
            break;
    }
}
*/
// No-op default implementation: return an empty vector
std::vector<double> ODEPar::SolveODE()
{
    return std::vector<double>(0);
}

bool ODEPar::Compare(const ODEPar rhs)
{
    if (numPars != rhs.numPars)
    {
        return false;
    }
    int sum = 0;
    for (size_t i = 0; i < numPars; ++i)
    {
        sum += (_pars[i] == rhs._pars[i] ? 1 : 0);
    }
    return sum == numPars;
}

ODEPar::ODEPar(int pars) : numPars(pars), _pars(pars, 1.0) {}
/*
ODEPar::ODEPar(int numPar, int numTrait, std::vector<double> traits, std::vector<double> pars) : numPars(numPar), numTraits(numTrait)
{
    _pars.resize(pars.size());
    for (size_t i = 0; i < numPars; ++i)
    {
        _pars[i] = pars[i];
    }

    _solutionTraits.resize(traits.size());
    for (size_t j = 0; j < numTraits; ++j)
    {
        _solutionTraits[j] = traits[j]; 
    }

}
*/


// Get an ODEPar from a vector of ODEPars
double ODEPar::getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec, bool incrementCount)
{
    for (auto& ODE : vec)
    {
        ODEPar source = *(ODE.get());
        if (target == source)
        {
            if (incrementCount) 
                ++*ODE;

            return source._AUC;
        }
    }
    return 0;
} 

std::vector<double> ODEPar::getPars(bool returnAUC)
{
    int n = numPars;
    int startIndex = 0;

    // If we're returning the AUC, make the vector bigger
    if (returnAUC)
    {
        startIndex = 1;
    }
    std::vector<double> result(n + startIndex);

    for (int i = 0; i < n; ++i)
    {
        result[i + startIndex] = _pars[i];
    }

    // First entry is AUC
    if (returnAUC)
    {
        result[0] = _AUC;
    }

    return result;
}

double ODEPar::getParValue(int i)
{
    return _pars[i];
}

void ODEPar::setParValue(int i, double val)
{
    _pars[i] = val;
    return;
}

void ODEPar::setParValue(std::vector<double> vals, bool firstElementIsAUC)
{
    int startIndex = 0;

    if (firstElementIsAUC)
    {
        startIndex++;
        setAUC(vals[0]);
    }    

    for (int i = startIndex; i < numPars; ++i)
    {
        setParValue(i, vals[i]);
    }
    return;
}

// Calculates area under the curve via the trapezoid method
#pragma omp declare simd
double ODEPar::AUC(const double &h, const double &a, const double &b)
{
	return ((a + b) * 0.5) * h;
}

// Calculates the response time [0], steady state concentration [1], and time to steady state [2] for a given ODE solution
std::vector<double> ODEPar::CalcSteadyState(const asc::Recorder &solution, const double& startTime, const double& stopTime, const int &solutionIndex)
{
    // First get the steady state and halfway point
    std::vector<double> result{0.0, 0.0, 0.0};
    double half = 0.0;
    static float epsilon = 0.001f;
    int steadyCount = 0;
    static int maxSteadyCount = 3;

    // Find start index in solution history
    int startIndex = (int)startTime * 10; // TODO: HACK: multiply by 10 because the sampling rate is 0.1
    int stopIndex = (int)stopTime * 10;

    // Make sure start index isn't out of bounds
    if (startIndex >= solution.history.size() - 2) {
        startIndex = solution.history.size() - 2;
    }

    if (stopIndex >= solution.history.size() - 1) {
        stopIndex = solution.history.size() - 1;
    }


    for (int i = startIndex + 1; i < stopIndex; ++i)
    {
        // concentrations for current and previous time points
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];
        if (std::abs(c2 - c1) < epsilon) {
            steadyCount++;
            if (steadyCount >= maxSteadyCount) {
                result[1] = c2; // steady state value
                result[2] = solution.history[i][0];
                break;
            }
            continue;
        }
        // Reset steadyCount if we're larger than the epsilon
        steadyCount = 0;
    }

    // Check that the result is valid
    if (std::isnan(result[2]) || result[2] > MAX_TIME) 
    {
        result[2] = 0.0;
    }

    if (std::isnan(result[1]) || result[1] > MAX_EXP) 
    {
        result[1] = MAX_EXP; // HACK: Big number, but not so big we run into floating point issues
    }

    // If there's no steady state, return early
    if (result[1] <= epsilon) {
        return result;
    }
    // Find the response time
    half = result[1] * 0.5;
    
    // Figure out where the halfway point is
    for (int i = startIndex + 1; i < solution.history.size() - 1; ++i)
    {
        double t1 = solution.history[i-1][0];
        double t2 = solution.history[i][0];
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];
        if ((c1 < half && c2 >= half) || (c1 > half && c2 <= half)) {
            result[0] = Interpolate(t1, c1, t2, c2, half) - startTime; // response time, relative to start time
            break;
        }
    }

    if (std::isnan(result[0]) || result[0] > MAX_TIME) 
    {
        result[0] = 0.0;
    }

    return result;
}

// Calculates response time to steady state [0], second steady state [1], time to steady state [2]
std::vector<double> ODEPar::CalcSecondSteadyState(const asc::Recorder &solution, const double& prevSteadyState, const double& prevSteadyStateTime, const int &solutionIndex)
{
    std::vector<double> result{0.0, 0.0, 0.0};
    double half = 0.0;
    static float epsilon = 0.001f;
    int steadyCount = 0;
    static int maxSteadyCount = 5;

    // Find start index in solution history
    int startIndex = (int)prevSteadyStateTime * 10; // TODO: HACK: multiply by 10 because the sampling rate is 0.1

    if (startIndex > solution.history.size() - 1) {
        startIndex = solution.history.size() - 1;
    }


    for (int i = startIndex + 1; i < solution.history.size() - 1; ++i)
    {
        // concentrations for current and previous time points
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];
        if (std::abs(c2 - c1) < epsilon) {
            steadyCount++;
            if (steadyCount >= maxSteadyCount) {
                result[1] = c2; // steady state value
                result[2] = solution.history[i][0];
                break;
            }
            continue;
        }
        // Reset steadyCount if we're larger than the epsilon
        steadyCount = 0;
    }

        // Check that the result is valid
    if (std::isnan(result[2]) || std::isinf(result[2])) 
    {
        result[2] = 0.0;
    }

    if (std::isnan(result[1]) || result[1] > MAX_EXP) 
    {
        result[1] = 10000; // HACK: Big number, but not so big we run into floating point issues
    }


    // Find the response time: taken from difference between old and new steady state
    half = std::abs(prevSteadyState - result[1]) * 0.5;

    // Figure out where the halfway point is
    for (int i = startIndex; i < solution.history.size(); ++i)
    {
        double t1 = solution.history[i-1][0];
        double t2 = solution.history[i][0];
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];
        if ((c1 < half && c2 >= half) || (c1 > half && c2 <= half)) {
            result[0] = Interpolate(t1, c1, t2, c2, half) - prevSteadyStateTime; // response time, relative to previous steady state
            break;
        }
    }

    if (std::isnan(result[0]) || std::isinf(result[0])) 
    {
        result[0] = 0.0;
    }

    return result;
}

// Finds the delay from a certain point before n ODE starts changing
// Since the PAR requires input from environment, we are looking for a change from the steady increase in Z
// afforded from this input (i.e. activation of the PAR/FFL circuit)
double ODEPar::CalcDelayTime(const asc::Recorder &solution, const double &startTime, const double &stopTime, const int &solutionIndex, const double &aZ, const double &baseline)
{
    double result = 0.0;

    float threshold = 0.0;

    int startIndex = (int)startTime * 10 + 1; // TODO: HACK: multiply by 10 because the sampling rate is 0.1

    int stopIndex = (int)stopTime * 10;

    // Find the threshold value: maximum change during simple regulation
    for (int i = 0; i < startIndex + 1; ++i) 
    {   
        float t1 = i/10.0;
        float t2 = (i+1)/10.0;
        // HACK: 0.1 is the sampling rate of the model
        float c1 = SimpleRegulation(t1, aZ, baseline);
        float c2 = SimpleRegulation(t2, aZ, baseline);
        float diff = std::abs(c2 - c1);

        if (diff > threshold)
        {
            threshold = diff;
        }
    }

    float epsilon = 0.001f;

    if (threshold < epsilon)
    {
        threshold = epsilon;
    }



    if (startIndex >= solution.history.size() - 1) {
        startIndex = solution.history.size() - 1;
    }

    if (stopIndex > solution.history.size() - 1) {
        stopIndex = solution.history.size() - 1;
    }

    for (int i = startIndex + 1; i < stopIndex; ++i)
    {
        double t = solution.history[i][0];
        double c1 = solution.history[i - 1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];

        double newDiff = std::abs(c2 - c1);

        // If the difference in concentration is greater than the epsilon, we've started moving and the delay is over
        if (newDiff > threshold) {
            result = t - startTime; // Offset by startTime
            break;
        }
    }

    return result;
}

// Maximum expression [0] and the response time to max expression (half maximum) [1]
std::vector<double> ODEPar::CalcMaxExpression(const asc::Recorder &solution, const float &startTime, const int &solutionIndex)
{
    std::vector<double> result {0.0, 0.0};
    double curMax = 0.0;
    double curTime = 0.0;

    int startIndex = (int)startTime * 10 + 1; // TODO: HACK: multiply by 10 because the sampling rate is 0.1

    if (startIndex > solution.history.size() - 1) {
        startIndex = solution.history.size() - 1;
    }

    for (int i = startIndex; i < solution.history.size(); ++i)
    {
        double curVal = solution.history[i][solutionIndex];
        if (curVal > curMax + 0.0001) {
            curMax = curVal;
            //curTime = solution.history[i][0];
        }
    }


    // Check values
    if (std::isnan(curMax) || curMax > MAX_EXP) 
    {
        curMax = 10000; // HACK: Big number, but not so big we run into floating point issues
    }

    float half = curMax * 0.5;

    // Get response time/half maximum expression
    for (int i = startIndex + 1; i < solution.history.size(); ++i)
    {
        double t1 = solution.history[i-1][0];
        double t2 = solution.history[i][0];
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];
        if ((c1 < half && c2 >= half) || (c1 > half && c2 <= half)) {
            curTime = Interpolate(t1, c1, t2, c2, half) - startTime; // response time, relative to startTime
            break;
        }
    }

    if (std::isnan(curTime) || curTime > MAX_TIME) 
    {
        curTime = 0.0;
    }

    result[0] = curMax;
    result[1] = curTime;

    return result;
}

// Amount of time spent above a threshold
double ODEPar::CalcTimeAboveThreshold(const asc::Recorder &solution, const double &threshold, const int &solutionIndex)
{
    static double DELTA = 0.1; // delta between two measurements
    double timeAboveThreshold = 0.0;

    for (int i = 0; i < solution.history.size(); ++i)
    {
        double curVal = solution.history[i][solutionIndex];

        if (curVal >= threshold) {
            timeAboveThreshold += DELTA;
        }
    }

    return timeAboveThreshold;
}