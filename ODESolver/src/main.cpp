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
        { "opt_file", required_argument,    0,  'O' },
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
    "-p             Specify the fitness optimum, the vector of trait values where fitness is\n"
    "               maximised. Separate traits with a comma ','.\n"
    "               Example: -p '2,3,5'\n"
    "-w             Specify the width of the multivariate Gaussian fitness function: larger\n"
    "               numbers mean stronger selection. Each value corresponds to the width on\n"
    "               a single trait axis. Separate values commas ','.\n"
    "               Example: -w 0.05,0.001,0.2\n"
    "-I             Specifies that the first column is an identifier for the combo.\n"
    "-O             Specify an input file that contains the optimum trait values and selection\n"
    "               strengths.\n"
    "               Example: -O '/path/to/optimum_width.csv'.\n"
    "               The first n columns should be the optimum values, where n is the number of\n"
    "               traits. The final columns should be the selection strength on each trait. If\n"
    "               -I is enabled, the first columns is the identifier row.\n"
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
        { "opt_file",   required_argument,  0,  'O' },
        {0,0,0,0}
    };

    int opt_idx = 0;
    int options;
    rapidcsv::Document doc;
    std::string output_path;
    std::string system_string;
    unsigned int nthreads = 1;
    std::string opt_str;

    // TODO: To start, we will only handle diagonal Sigma matrices, so we can feed in a vector
    // of variances instead of a full v/cov matrix
    std::string width_str;
    int id = 0;

    // Handle input optima
    rapidcsv::Document doc_opt;

    while (options != -1)
    {
        options = getopt_long(argc, argv, "hi:o:s:t:p:w:IO:", voptions, &opt_idx);
    
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
            nthreads = (std::stoi(optarg) > 0) ? std::stoi(optarg) : 1;
            continue;

        case 'p':
            opt_str = (std::string)optarg;
            continue;
        
        case 'w':
            width_str = (std::string)optarg;
            continue;

        case 'I':
            id = 1;
            continue;
        
        case 'O':
            doc_opt.Load(((std::string)optarg), rapidcsv::LabelParams(-1, -1));
            continue;
        
        case -1:
            break;
        }
    
    }

    // Handle width and optimum conversion from string input to std::vector<double>
    std::vector<double> optimum;
    std::vector<double> width;

    if (doc_opt.GetRowCount() == 0)
    {
        std::string temp_str;
        std::stringstream ss(opt_str);
        
        if (opt_str == "") 
        {
            optimum.emplace_back(0.0);
        }
        else 
        {
            // Fill the optimum vector
            while (std::getline(ss, temp_str, ',')) 
            {
                // Remove any spaces
                temp_str.erase(std::remove_if(temp_str.begin(), temp_str.end(), 
                [](const char& c) {return &c == " ";} ), temp_str.end());
                
                // Add to end
                optimum.emplace_back(std::stod(temp_str));
            }
            
        }
        
        if (width_str == "")
        {
            width.emplace_back(1.0);
        }
        else
        {
            // Fill the width vector
            ss = std::stringstream(width_str);
            while (std::getline(ss, temp_str, ',')) 
            {
                // Remove any spaces
                temp_str.erase(std::remove_if(temp_str.begin(), temp_str.end(), 
                [](const char& c) {return &c == " ";} ), temp_str.end());
                
                width.emplace_back(std::stod(temp_str));
            }
        }
        
        // If the width and optimum vectors are different sizes, something has gone wrong
        if (width.size() != optimum.size()) 
        {
            std::cerr << "Optimum and width are unequal sizes!" << std::endl;
            return 1;
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
    const auto SolveODESystem = [&doc, &result, &width, &optimum, &doc_opt, &id, &motif](const unsigned i)
    {
        // Get molecular trait values
        const std::vector<double> parameters = doc.GetRow<double>(i);

        // Get id - may not be set, but we'll check later
        const double id_num = parameters[0];

        // Solve for phenotype: offset by 1 if id is true
        std::unique_ptr ODE = ODEPar::MakeODEPtr(motif);

        // Check if the row has the right number of elements
        if (parameters.size() - id != ODE->GetNumPars())
        {
            std::cerr << "Error in parsing row " << i << 
            ": number of input parameters does not match ODE parameter count." << std::endl;
            return;
        }

        // If we have an id, we need to offset parameters by 1 since first entry is id
        // Set ODE parameters
        //std::cout << "Setting par values with id = " << id << std::endl;
        ODE->setParValue({parameters.begin() + id, parameters.end()}, false);

        //std::cout << "Solving ODE" << std::endl;

        // Solve ODE
        ODE->SolveODE();
        
        std::vector<double> thisRunOpt;
        std::vector<double> thisRunWidth;

        // Figure out fitness for input data
        //std::cout << "Opt file: " << doc_opt.GetRowCount() << std::endl;
        //std::cout << "Input file: " << doc.GetRowCount() << std::endl;

        if (doc_opt.GetRowCount() == doc.GetRowCount()) {
            //std::cout << "Getting optima from file..." << std::endl;
            size_t n = ODE->GetNumTraits();
            //std::cout << n << std::endl;
            thisRunOpt.reserve(n);
            thisRunWidth.reserve(n);
            const std::vector<double> opt_width = doc_opt.GetRow<double>(i);
            // Get indexes for each 
            //std::cout << "Traits = " << n << std::endl;
            auto startOptIndex = opt_width.begin() + id;
            auto endOptIndex = startOptIndex + n;
            auto endWidthIndex = opt_width.end();
            thisRunOpt.insert(thisRunOpt.begin(), startOptIndex, endOptIndex);
            thisRunWidth.insert(thisRunWidth.begin(), endOptIndex, endWidthIndex);
        }

        if (doc_opt.GetRowCount() == 0) {
            thisRunOpt = optimum;
            thisRunWidth = width;
        }

        //std::cout << "Optimum: " << thisRunOpt.size() << std::endl;
        //std::cout << "Width: " << thisRunWidth.size() << std::endl;

        // Write to output vector
        result[i] = std::make_unique<std::string>(ODE->printPars(thisRunWidth, thisRunOpt, (char* const)","));
        // if bool, we need to attach parameter[0] to the front
        if (id > 0)
        {
            result[i].get()->insert(0, std::to_string(id_num) + ",");
        }

        return;
    };
    
    // Parallel compute for each row in input
    //SolveODESystem(0);
    
    //////////////
    parallelutil::queue_based_parallel_for(doc.GetRowCount(), SolveODESystem, 1);//nthreads);

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