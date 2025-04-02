// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include "derivatives.h"
#include "schemes.h"
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


JsonData json_reader(const std::string& filename) {
    JsonData config;
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

    return config;
}


template <int Precision>
void runSolver(JsonData config) {
    const int domain_size = config.integerValues["domain_size"];

    using Real = typename PrecisionToType<Precision>::Type;

    std::cout << "--- Using " << (Precision == 32 ? "float" : (Precision == 64 ? "double" : "long double")) << " (Second Derivatives) ---" << std::endl;
    Real dx = static_cast<Real>(0.1);
    Real x0 = static_cast<Real>(1.0);
    Real f_il3 = std::sin(x0 - 3 * dx);
    Real f_il2 = std::sin(x0 - 2 * dx);
    Real f_il1 = std::sin(x0 - dx);
    Real f_i0 = std::sin(x0);
    Real f_ir1 = std::sin(x0 + dx);
    Real f_ir2 = std::sin(x0 + 2 * dx);
    Real f_ir3 = std::sin(x0 + 3 * dx);

    std::cout << "Exact second derivative of sin(x) at x=" << x0 << " is -sin(" << x0 << ") = " << -std::sin(x0) << std::endl;
    std::cout << "FDS1 (O(h)): " << fds1_second_derivative(dx, f_i0, f_ir1, f_ir2) << std::endl;
    std::cout << "BDS1 (O(h)): " << bds1_second_derivative(dx, f_il2, f_il1, f_i0) << std::endl;
    std::cout << "CDS2 (O(h^2)): " << cds2_second_derivative(dx, f_il1, f_i0, f_ir1) << std::endl;
    std::cout << "FDS2 (O(h^2)): " << fds2_second_derivative(dx, f_i0, f_ir1, f_ir2, f_ir3) << std::endl;
    std::cout << "BDS2 (O(h^2)): " << bds2_second_derivative(dx, f_il3, f_il2, f_il1, f_i0) << std::endl;
    std::cout << std::endl;


    //-------------------------------------------------------------------------
    // Raw pointer version
    Real* input = create_initialized_array<Real, Real>(domain_size, 1.0);

    // Compute derivatives
    Real* derivatives = central_diff_first_derivative(input, domain_size, dx);

    // Use the derivatives...
    for(int i = 0; i < domain_size; ++i) {
        std::cout << "Value at index " << i << ": " << input[i] << '\n';
    }


    // Use the derivatives...
    for(int i = 0; i < domain_size; ++i) {
        std::cout << "Derivative at index " << i << ": " << derivatives[i] << '\n';
    }

    // Cleanup memory
    delete[] derivatives;
    //-------------------------------------------------------------------------
}



int main() {
    JsonData config;
    std::string filename = "configuration.json"; // Replace with your JSON filename

    config = json_reader(filename);

    try {
        simple_demo();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }


    int precision = config.integerValues["float_precision"];

    if (precision == 32) {
        runSolver<32>(config);
    } else if (precision == 64) {
        runSolver<64>(config);
    } else if (precision == 128) {
        runSolver<128>(config);
    } else {
        std::cout << "Error: Unsupported float precision: " << precision << std::endl;
        return 1;
    }





    



    return 0;
}