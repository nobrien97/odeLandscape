#include "ascent/Ascent.h"
#include "include/rapidcsv.h"
#include "ODESystem.h"
#include "odePar.h"
#include "getopt.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "include/parallel-util.hpp"

#define no_arg 0
#define arg_required 1
#define arg_optional 2

    const struct option longopts[] =
    {
        { "input",    required_argument,        0,  'i' },
        { "output",   required_argument,        0,  'o' },
        { "help",           no_argument,        0,  'h' },
        { "threads",  required_argument,        0,  't' }, 
        {0,0,0,0}
    };


void doHelp(char* appname) {
    std::fprintf(stdout,
    "ODE Landscaper: A program for calculating NAR ODE parameter combinations.\n"
    "\n"
    "This program generates phenotypes from a list of molecular trait values.\n"
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
    "-t             Specify number of threads to use while calculating values.\n"
    "               Example: -t 4\n"
    "\n",
    appname,
    appname
    );
}


int main(int argc, char* argv[])
{
    const struct option voptions[] = 
    {
        { "help",        no_argument,        0,  'h' },
        { "input",       required_argument,  0,  'i' },
        { "output",      required_argument,  0,  'o' },
        { "threads",     required_argument,  0,  't' },
        {0,0,0,0}
    };

    int opt_idx = 0;
    int options;
    rapidcsv::Document doc;
    size_t length = 1000;
    std::string output_path;
    unsigned nthreads = 1;

    while (options != -1)
    {
        options = getopt_long(argc, argv, "hi:o:t:", voptions, &opt_idx);
    
        switch (options)
        {
        case 'h':
            doHelp(argv[0]);
            return 0;

        case 'l':
            length = (size_t)optarg;
            continue;

        case 'i':
            //TODO: Read csv file, fill parameter ranges
            doc.Load(((std::string)optarg));
            continue;
        
        case 'o':
            //TODO: write output .csv with parameter combinations and phenotypes
            output_path = (std::string)optarg;
            continue;

        case 't':
            nthreads = (unsigned)optarg;
            continue;
        
        case -1:
            break;
        }
    
    }

    // Check nthreads is <= number cores
    const auto proc_count = std::thread::hardware_concurrency();
    nthreads = (nthreads <= proc_count) ? nthreads : proc_count;

// Lambda to solve NAR system
    const auto SolveNAR = [&doc](const unsigned i)
    {
        // Get molecular trait values
        const std::vector<double> parameters = doc.GetRow<double>(i);

        // Solve for phenotype
        ODESystem ODE(ODEPar(-1.0, parameters[0], parameters[1], parameters[2], parameters[3]));
        ODE.calculatePhenotype();

        // Write to a file on this thread
        std::ofstream cur_file;
        cur_file.open("solution_" + std::to_string(i)  + ".csv");
        cur_file << ODE.printPars(",") + "\n";
        cur_file.close();

        return;
    };
    
    // Setup output vector
    std::vector<std::unique_ptr<std::string>> result(doc.GetRowCount(), nullptr);

    

    // Parallel compute for each row in input
    parallelutil::queue_based_parallel_for(doc.GetRowCount(), SolveNAR, nthreads);


    /* TODO: thread here - each thread:
        Gets a parameter combination
        Runs the parameter combination
        Writes to file using some synchronized file writing class to avoid race conditions 
    */ 


    return 0;
}