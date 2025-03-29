// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include <iostream>

using namespace Cantera;

// The actual code is put into a function that can be called from the main program.
void simple_demo()
{
    // Create a new Solution object
    auto sol = newSolution("h2o2.yaml");
    auto gas = sol->thermo();

    // Set the thermodynamic state by specifying T (500 K) P (2 atm) and the mole
    // fractions. Note that the mole fractions do not need to sum to 1.0 - they will
    // be normalized internally. Also, the values for any unspecified species will be
    // set to zero.
    gas->setState_TPX(500.0, 2.0*OneAtm, "H2O:1.0, H2:8.0, AR:1.0");

    // Print a summary report of the state of the gas.
    std::cout << gas->report() << std::endl;
}


int main() {
    JsonData config;
    std::string filename = "configuration.json"; // Replace with your JSON filename

    if (readAndParseJson(filename, config)) {
        std::cout << "configuration.json file parsed successfully!" << std::endl;

        // Access the parsed values:
        if (config.stringValues.count("data_folder01")) {
            std::cout << "Path to cantera files: " << config.stringValues["data_folder01"] << std::endl;
        }

        if (config.integerValues.count("run_num")) {
            std::cout << "Run number chosen: " << config.integerValues["run_num"] << std::endl;
        }

        if (config.integerValues.count("float_precision")) {
            std::cout << "Precision chosen for the run: " << config.integerValues["float_precision"] << std::endl;
        }

        if (config.stringValues.count("anotherString")) {
            std::cout << "Another String Value: " << config.stringValues["anotherString"] << std::endl;
        }

        if (config.floatValues.count("floatValue")) {
            std::cout << "Float Value: " << config.floatValues["floatValue"] << std::endl;
        }

        // Access other values similarly based on your JSON structure.
    } else {
        std::cerr << "Failed to parse JSON file." << std::endl;
    }


    try {
        simple_demo();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }


    return 0;
}