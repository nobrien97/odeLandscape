#include "ascent/Ascent.h"
#include "include/rapidcsv.h"
#include "ODESystem.h"
#include "odePar.h"
#include "getopt.h"
#include "matrixClass.h"
#include <vector>
#include <iostream>

#define no_arg 0
#define arg_required 1
#define arg_optional 2

    const struct option longopts[] =
    {
        { "input",    required_argument,        0,  'i' },
        { "output",   required_argument,        0,  'o' },
        { "help",           no_argument,        0,  'h' },
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
    "               Input should contain two columns titled \'min\' and \'max'\n"
    "               with each row specifying a minimum and maximum value for a\n"
    "               molecular trait. The program will calculate values across\n"
    "               entire range. For a single molecular trait combination,\n"
    "               specify a single column named \'value\'.\n"
    "\n"
    "-o             Specify output filepath.\n"
    "               Example -o /path/to/output.csv\n"
    "\n",
    appname,
    appname
    );
}

std::vector<std::vector<double>> getDocValues(const rapidcsv::Document& doc, size_t seq_length)
{
    size_t col_count = doc.GetColumnCount();
    size_t row_count = doc.GetRowCount();
    
    if (col_count == 1)
    {
        std::vector<std::vector<double>> result = std::vector<std::vector<double>>(row_count, std::vector<double>(1));
        
        for (int i = 0; i < row_count; ++i)
        {
            result[i][0] = doc.GetCell<double>(0, i+1);
        }
        return result;
    }

    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(row_count);
    
    for (int i = 0; i < row_count; ++i)
    {
        double min = doc.GetCell<double>(0, i+1);
        double max = doc.GetCell<double>(1, i+1);
        double by = (max - min)/seq_length;

        std::vector<double> row_values = seq(min, max, by);
        
        result[i] = row_values;
    }

    return result;
}

std::vector<double> seq(double from, double to, double by)
{
    int length = (int)((to - from)/by);
    std::vector<double> result;
    for (int i = 0; i < length; ++i)
    {
        result.emplace_back(from + i*by);
    }
    return result;
}

int main(int argc, char* argv[])
{
    const struct option voptions[] = 
    {
        { "help",        no_argument,        0,  'h' },
        { "input",       required_argument,  0,  'i' },
        { "output",      required_argument,  0,  'o' },
        {0,0,0,0}
    };

    int opt_idx = 0;
    int options;
    rapidcsv::Document doc;
    size_t length = 1000;
    std::string output_path;
    std::vector<std::vector<double>> sequences;

    while (options != -1)
    {
        options = getopt_long(argc, argv, "hi:o:", voptions, &opt_idx);
    
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
            sequences = getDocValues(doc, length);
            continue;
        
        case 'o':
            //TODO: write output .csv with parameter combinations and phenotypes
            output_path = optarg;
            continue;
        
        case -1:
            break;
        }
    
    }

    // TODO: Fix this, get command line input working - make it work on a file and extract the values from there, 
    // operate over a range (e.g. specify aZ min to aZ max)
    std::vector<double> aZ = sequences[0];

    /* TODO: thread here - each thread:
        Gets a parameter combination
        Runs the parameter combination
        Writes to file using some synchronized file writing class to avoid race conditions 
    */ 

    ODEPar parameters(-1.0, );

    ODESystem ODE(parameters);

    // Calculate phenotype, write to file
    ODE.calculatePhenotype();

    return 0;
}