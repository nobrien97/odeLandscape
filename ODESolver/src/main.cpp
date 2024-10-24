#include "ascent/Ascent.h"
#include "version.h"
#include "rapidcsv.h"
#include "odePar.h"
#include "getopt.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "parallel-util.hpp"

#define no_arg 0
#define arg_required 1
#define arg_optional 2

    const struct option longopts[] =
    {
        { "input",    required_argument,    0,  'i' },
        { "output",   required_argument,    0,  'o' },
        { "system",   required_argument,    0,  's' },
        { "help",           no_argument,    0,  'h' },
        { "threads",  required_argument,    0,  't' }, 
        { "optimum",  required_argument,    0,  'p' },
        { "width",    required_argument,    0,  'w' },
        {0,0,0,0}
    };


void doHelp(char* appname) {
    std::fprintf(stdout,
    "ODE Landscaper (V%i.%i): A program for calculating NAR ODE parameter combinations.\n"
    "\n"
    "This program generates phenotype and fitness landscapes from a list of molecular\n"
    "trait values.\n"
    "Usage: %s [OPTION]...\n"
    "Example: %s -h\n"
    "\n"
    "-h             Print this help manual.\n"
    "\n"
    "-i             Specify input filepath.\n"
    "               Example: -i /path/to/file.csv\n"
    "               Input should contain a column for each molecular trait\n"
    "               with each row specifying a parameter combination.\n"
    "\n"
    "-o             Specify output filepath.\n"
    "               Example: -o /path/to/output.csv\n"
    "\n"
    "-s             Specify the ODE system to use, PAR, NAR, FFLC1, FFLI1, or FFBH.\n"
    "               Example: -s 'PAR'\n"
    "\n"
    "-t             Specify number of threads to use while calculating values.\n"
    "               Example: -t 4\n"
    "\n"
    "-p             Specify the fitness optimum, the phenotype where fitness is\n"
    "               maximised.\n"
    "               Example: -p 2\n"
    "-w             Specify the width of the Gaussian fitness function: larger\n"
    "               numbers mean stronger selection.\n"
    "               Example: -w 0.05\n"
    "-I             Specifies that the first column is an identifier for the combo.\n"
    "\n",
    ODELandscaper_VERSION_MAJOR,
    ODELandscaper_VERSION_MINOR,
    appname,
    appname
    );
}


int main(int argc, char* argv[])
{
    const struct option voptions[] = 
    {
        { "help",       no_argument,        0,  'h' },
        { "input",      required_argument,  0,  'i' },
        { "output",     required_argument,  0,  'o' },
        { "system",     required_argument,  0,  's' },
        { "threads",    required_argument,  0,  't' },
        { "optimum",    required_argument,  0,  'p' },
        { "width",      required_argument,  0,  'w' },
        { "identifier", no_argument,        0,  'I' },
        {0,0,0,0}
    };

    int opt_idx = 0;
    int options;
    rapidcsv::Document doc;
    std::string output_path;
    std::string system_string;
    unsigned int nthreads = 1;
    double optimum = 0;
    double width = 1;
    int id = 0;

    while (options != -1)
    {
        options = getopt_long(argc, argv, "hi:o:s:t:p:w:I", voptions, &opt_idx);
    
        switch (options)
        {
        case 'h':
            doHelp(argv[0]);
            return 0;

        case 'i':
            // Read csv file
            doc.Load(((std::string)optarg), rapidcsv::LabelParams(-1, -1));
            continue;
        
        case 'o':
            // write output .csv with parameter combinations and phenotypes
            output_path = (std::string)optarg;
            continue;

        case 's':
            // choose ODE system
            system_string = (std::string)optarg;
            continue;

        case 't':
            std::cout << "threads: " << optarg << std::endl;
            nthreads = (std::stoi(optarg) > 0) ? std::stoi(optarg) : 1;
            continue;

        case 'p':
            optimum = std::stod(optarg);
            continue;
        
        case 'w':
            width = std::stod(optarg);
            continue;

        case 'I':
            id = 1;
            continue;
            
        case -1:
            break;
        }
    
    }

    if (doc.GetRowCount() == 0) 
    {
        std::cerr << "Please provide a valid input .csv file with \"-i path/to/file.csv\"." << std::endl;
        return 1;
    }

    // Check nthreads is <= number cores
    const auto proc_count = std::thread::hardware_concurrency();
    nthreads = (nthreads <= proc_count) ? nthreads : proc_count;

    // Setup output vector
    std::vector<std::unique_ptr<std::string>> result(doc.GetRowCount());

    // Get motif type
    ODEPar::motif_enum motif = ODEPar::HashMotifString(system_string);

// Lambda to solve NAR system
    const auto SolveODESystem = [&doc, &result, &width, &optimum, &id, &motif](const unsigned i)
    {
        // Get molecular trait values
        const std::vector<double> parameters = doc.GetRow<double>(i);

        // Get id - may not be set, but we'll check later
        const double id_num = parameters[0];

        // Solve for phenotype: offset by 1 if id is true
        std::unique_ptr ODE = ODEPar::MakeODEPtr(motif);

        // Check if the row has the right number of elements
        if (parameters.size() - id != ODE->numPars)
        {
            std::cout << "Error in parsing row " << i << 
            ": number of input parameters does not match ODE parameter count." << std::endl;
            return;
        }

        // If we have an id, we need to offset parameters by 1 since first entry is id
        // Set ODE parameters
        ODE->setParValue({parameters.begin() + id, parameters.end()}, false);

        // Solve ODE
        ODE->SolveODE();

        // Write to output vector
        result[i] = std::make_unique<std::string>(ODE->printPars(width, optimum, (char* const)","));
        // if bool, we need to attach parameter[0] to the front
        if ((bool)id)
        {
            result[i].get()->insert(0, std::to_string(id_num) + ",");
        }

        return;
    };
    
    // Parallel compute for each row in input
    parallelutil::queue_based_parallel_for(doc.GetRowCount(), SolveODESystem, nthreads);

    // Write to file
    std::ofstream file;
    file.open(output_path);

    for (int i = 0; i < result.size(); ++i)
    {
        file << (*result[i]) << "\n";
    }

    file.close();

    return 0;
}