// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include "derivatives.h"
#include "schemes.h"
#include "FluidSolver1D.h"
#include <iostream>

using namespace Cantera;


template <typename T>
std::shared_ptr<Cantera::ThermoPhase> gas_h2o2_creator(T temp, T pressure, std::string mole_fractions){
    
    std::shared_ptr<Cantera::Solution> sol = newSolution("h2o2.yaml");
    std::shared_ptr<Cantera::ThermoPhase> gas_h2o2 = sol->thermo();

    gas_h2o2->setState_TPX(temp, pressure, mole_fractions);

    return gas_h2o2;
}



JsonData json_reader(const std::string& filename) {
    JsonData config;
    if (readAndParseJson(filename, config)) {
        std::cout << "configuration.json file parsed successfully!" << std::endl;
    } else {
        std::cerr << "Failed to parse JSON file." << std::endl;
    }
    return config;
}



template <int Precision>
void FluidSolver1DFunc(JsonData config) {
    using Real = typename PrecisionToType<Precision>::Type;
    
    const std::string data_folder01 = config.stringValues["data_folder01"];
    const int run_num = config.integerValues["run_num"];
    const int float_precision = config.integerValues["float_precision"];
    const int boundary_mode = config.integerValues["boundary_mode"];
    const double init_temp = config.floatValues["init_temp"];
    const double init_pressure = config.floatValues["init_pres"];
    const int Nx = config.integerValues["Nx"];
    const double domain_length = config.floatValues["domain_length"];
    const double timestep = 5e-9; //static_cast<double>(config.floatValues["timestep"]);
    const int n_iters_total = config.integerValues["n_iters_total"];
    const int n_iters_save = config.integerValues["n_iters_save"];
    const std::string scheme_to_use = config.stringValues["scheme_to_use"];

    Real T0 = static_cast<Real>(init_temp);
    Real P0 = static_cast<Real>(init_pressure);

    std::shared_ptr<Cantera::ThermoPhase> gas_h2o2 = gas_h2o2_creator(T0, P0, air_gas_comp);

    FluidSolver1D solver = FluidSolver1D<Real>(gas_h2o2, T0, P0, 
        fuel_gas_comp, static_cast<Real>(domain_length), Nx, static_cast<Real>(timestep),
        scheme_to_use);
    
    std::vector<Real> rho_new;
    std::vector<Real> rho_u_new;
    std::vector<Real> rho_e0_new;
    std::vector<std::vector<Real>> rho_ys_new;

    int case1 = 1;

    std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new) = solver.get_initial_vals(case1);
    
    std::cout << "so far, so good" << std::endl;
}


int main() {
    JsonData config;
    std::string filename = "configuration.json"; // Replace with your JSON filename

    config = json_reader(filename);
    int precision = config.integerValues["float_precision"];

    try{
        FluidSolver1DFunc<64>(config);
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }

    return 0;
}