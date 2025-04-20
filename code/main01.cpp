// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include "log_heap.h"
#include "derivatives.h"
#include "FluidSolver1D.h"
#include <iostream>
#include <chrono>
#include <stdexcept>

using namespace Cantera;


template <typename T>
std::shared_ptr<Cantera::Solution> gas_h2o2_creator(T temp, T pressure, std::string mole_fractions){
    std::shared_ptr<Cantera::Solution> gas_h2o2 = newSolution("../data/cantera_files/h2o2.yaml", "ohmech");
    gas_h2o2->thermo()->setState_TPX(temp, pressure, mole_fractions);

    return gas_h2o2;
}


JsonData json_reader(const std::string& filename) {
    JsonData config;
    if (readAndParseJson(filename, config)) {
        std::cout << "configuration.json file parsed successfully!" << std::endl;
        std::cout << "*********************************************************" << std::endl;
        std::cout << "The options chosen for the run are:" << std::endl;
        std::cout << "*********************************************************" << std::endl;
        std::cout << "data_folder01:      " << config.stringValues["data_folder01"] << std::endl;
        std::cout << "run_num:            " << config.stringValues["run_num"] << std::endl;
        std::cout << "float_precision:    " << config.integerValues["float_precision"] << std::endl;
        std::cout << "boundary_mode:      " << config.integerValues["boundary_mode"] << std::endl;
        std::cout << "init_temp:          " << config.floatValues["init_temp"] << std::endl;
        std::cout << "init_pres:          " << config.floatValues["init_pres"] << std::endl;
        std::cout << "Nx:                 " << config.integerValues["Nx"] << std::endl;
        std::cout << "domain_length:      " << config.floatValues["domain_length"] << std::endl;
        std::cout << "timestep:           " << config.floatValues["timestep"] << std::endl;
        std::cout << "n_iters_total:      " << config.integerValues["n_iters_total"] << std::endl;
        std::cout << "n_iters_save:       " << config.integerValues["n_iters_save"] << std::endl;
        std::cout << "scheme_to_use:      " << config.stringValues["scheme_to_use"] << std::endl;
        std::cout << "case # of problem:  " << config.integerValues["case_problem"] << std::endl;
        std::cout << "*********************************************************" << std::endl;

        if (config.integerValues["boundary_mode"] == 1) {
            std::cout << "Periodic boundary condition chosen " << std::endl;
        } else if (config.integerValues["boundary_mode"] == 2) {
            std::cout << "Both boundary condition chosen as Dirchlet" << std::endl;
            std::cout << "Left boundary value:  " << config.floatValues["left_boundary_dirichlet"] << std::endl;
            std::cout << "Right boundary value: " << config.floatValues["right_boundary_dirichlet"] << std::endl;
        } else if (config.integerValues["boundary_mode"] == 3) {
            std::cout << "Both boundary condition chosen as Neumann" << std::endl;
            std::cout << "Left boundary value:  " << config.floatValues["left_boundary_neumann"] << std::endl;
            std::cout << "Right boundary value: " << config.floatValues["right_boundary_neumann"] << std::endl;
        } else if (config.integerValues["boundary_mode"] == 4) {
            std::cout << "Left boundary condition -> Neumann, right boundary condition -> Dirichlet" << std::endl;
            std::cout << "Left boundary value:  " << config.floatValues["left_boundary_neumann"] << std::endl;
            std::cout << "Right boundary value: " << config.floatValues["right_boundary_dirichlet"] << std::endl;
        } else if (config.integerValues["boundary_mode"] == 5) {
            std::cout << "Left boundary condition -> Dirichlet, right boundary condition -> Neumann" << std::endl;
            std::cout << "Left boundary value:  " << config.floatValues["left_boundary_dirichlet"] << std::endl;
            std::cout << "Right boundary value: " << config.floatValues["right_boundary_neumann"] << std::endl;
        }

        std::cout << "*********************************************************" << std::endl;
    } else {
        std::cerr << "Failed to parse JSON file." << std::endl;
    }
    return config;
}



template <int Precision>
void FluidSolver1DFunc(JsonData config) {
    std::vector<std::chrono::steady_clock::time_point> time_checkpoints_start(0);
    std::vector<std::chrono::steady_clock::time_point> time_checkpoints_end(0);
    time_checkpoints_start.push_back(std::chrono::steady_clock::now());
    
    using Real = typename PrecisionToType<Precision>::Type;
    
    std::string data_folder01 = config.stringValues["data_folder01"];
    std::string outputs_folder01 = config.stringValues["outputs_folder01"];
    std::string run_num = config.stringValues["run_num"];
    const int float_precision = config.integerValues["float_precision"];

    // Boundary conditions
    const int boundary_mode = config.integerValues["boundary_mode"];
    Real left_boundary_dirichlet = static_cast<Real>(config.floatValues["left_boundary_dirichlet"]);
    Real right_boundary_dirichlet = static_cast<Real>(config.floatValues["right_boundary_dirichlet"]);
    Real left_boundary_neumann = static_cast<Real>(config.floatValues["left_boundary_neumann"]);
    Real right_boundary_neumann = static_cast<Real>(config.floatValues["right_boundary_neumann"]);

    // Initial conditions
    Real init_temp = static_cast<Real>(config.floatValues["init_temp"]);
    Real init_pressure = static_cast<Real>(config.floatValues["init_pres"]);

    // Domain conditions
    const int Nx = config.integerValues["Nx"];
    Real domain_length = static_cast<Real>(config.floatValues["domain_length"]);

    // Time parameters
    Real timestep = static_cast<Real>(static_cast<double>(config.floatValues["timestep"]));
    int n_iters_total = config.integerValues["n_iters_total"];
    int n_iters_save = config.integerValues["n_iters_save"];
    std::string scheme_to_use = config.stringValues["scheme_to_use"];
    int case1 = config.integerValues["case_problem"];

    Real T0 = static_cast<Real>(init_temp);
    Real P0 = static_cast<Real>(init_pressure);
    const std::string air_comp = gas_composition("air");
    const std::string fuel_comp = gas_composition("fuel");
    const std::string fuel_comp2 = gas_composition("fuel2");
    const std::string sample = gas_composition("other");


    // ******************************************************************
    // ******************************************************************
    // ******************************************************************

    std::shared_ptr<Cantera::Solution> gas_h2o2;
    if (case1 == 1){
        gas_h2o2 = gas_h2o2_creator(T0, P0, fuel_comp);
    } else if (case1 == 3){
        gas_h2o2 = gas_h2o2_creator(T0, P0, fuel_comp2);
    } else if (case1 == 2) {
        gas_h2o2 = gas_h2o2_creator(T0, P0, fuel_comp);
    } else {
        throw std::invalid_argument("Invalid case selected: " + std::to_string(case1) + ". Only case 1, 2, 3 available");
    }
    

    FluidSolver1D solver = FluidSolver1D<Real>(gas_h2o2, T0, P0, fuel_comp, domain_length, Nx, timestep,
                                                scheme_to_use, boundary_mode, left_boundary_dirichlet, right_boundary_dirichlet,
                                                left_boundary_neumann, right_boundary_neumann);
    
    std::vector<Real> rho_new, rho_u_new, rho_e0_new;
    std::vector<std::vector<Real>> rho_ys_new;
    std::vector<Real> temperature_old, pressure_old;

    std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new) = solver.get_initial_vals(case1);

    time_checkpoints_end.push_back(std::chrono::steady_clock::now());
    
    std::vector<std::string> variables_to_save = {"u", "dens", "temp", "pres", "h2", "h", "o", "o2", "oh", "h2o", "ho2", "h2o2", "ar", "n2"};
    std::vector<std::string> filenames(0);
    
    for (string variable : variables_to_save) {
        filenames.push_back(make_filename(outputs_folder01, variable, run_num));    
    }
    
    for (string filename : filenames) {
        clear_file(filename);
    }

    for (int i=0; i<n_iters_total; i++){
        time_checkpoints_start.push_back(std::chrono::steady_clock::now());
        
        Real time_keeper = static_cast<Real>(i) * timestep;

        std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old, pressure_old) = solver.solver_func(rho_new, rho_u_new, rho_e0_new, rho_ys_new);
        std::vector<Real> u_values = divide_vectors(rho_u_new, rho_new);
        std::vector<std::vector<Real>> ys_values = divide_vector_2nd_order_1st_order<Real>(rho_ys_new, rho_new);
        std::vector<Real> h2_vals = column_returner(ys_values, 0);
        std::vector<Real> h_vals = column_returner(ys_values, 1);
        std::vector<Real> o_vals = column_returner(ys_values, 2);
        std::vector<Real> o2_vals = column_returner(ys_values, 3);
        std::vector<Real> oh_vals = column_returner(ys_values, 4);
        std::vector<Real> h2o_vals = column_returner(ys_values, 5);
        std::vector<Real> ho2_vals = column_returner(ys_values, 6);
        std::vector<Real> h2o2_vals = column_returner(ys_values, 7);
        std::vector<Real> ar_vals = column_returner(ys_values, 8);
        std::vector<Real> n2_vals = column_returner(ys_values, 9);

        // print_vector(u_values);
        // print_vector(temperature_old);
        // std::cout << "Iteration: " << i << " running." << std::endl;

        time_checkpoints_end.push_back(std::chrono::steady_clock::now());

        if (i % n_iters_save == 0){
            file_writer_vector(time_keeper, u_values, filenames[0], 0);
            file_writer_vector(time_keeper, rho_new, filenames[1], 0);
            file_writer_vector(time_keeper, temperature_old, filenames[2], 0);
            file_writer_vector(time_keeper, pressure_old, filenames[3], 1);
            file_writer_vector(time_keeper, h2_vals, filenames[4], 0);
            file_writer_vector(time_keeper, h_vals, filenames[5], 0);
            file_writer_vector(time_keeper, o_vals, filenames[6], 0);
            file_writer_vector(time_keeper, o2_vals, filenames[7], 0);
            file_writer_vector(time_keeper, oh_vals, filenames[8], 0);
            file_writer_vector(time_keeper, h2o_vals, filenames[9], 0);
            file_writer_vector(time_keeper, ho2_vals, filenames[10], 0);
            file_writer_vector(time_keeper, h2o2_vals, filenames[11], 0);
            file_writer_vector(time_keeper, ar_vals, filenames[12], 0);
            file_writer_vector(time_keeper, n2_vals, filenames[13], 0);
        }
    }

    double total_time = 0.0;
    for (int i=0; i<time_checkpoints_start.size(); i++){
        // std::cout << "Loop time: " << std::chrono::duration_cast<std::chrono::microseconds>(time_checkpoints_end[i] - time_checkpoints_start[i]).count() /1000000.0 << std::endl;
        total_time += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(time_checkpoints_end[i] - time_checkpoints_start[i]).count() /1000000.0);
    }
    std::cout << "Running time (without I/O time) is: " << total_time << " seconds" << std::endl;

    std::cout << "so far, so good" << std::endl;
}


int main() {

    std::chrono::steady_clock::time_point begin_time_checkpoint1 = std::chrono::steady_clock::now();
    
    JsonData config;
    std::string filename = "configuration.json"; // Replace with your JSON filename

    config = json_reader(filename);
    int precision = config.integerValues["float_precision"];

    try{
        if (precision == 32){
            std::cout << "Selected precision 32" << std::endl;
            FluidSolver1DFunc<32>(config);
        } else if (precision == 64){
            std::cout << "Selected precision 64" << std::endl;
            FluidSolver1DFunc<64>(config);
        } else if (precision == 128){
            std::cout << "Selected precision 128" << std::endl;
            FluidSolver1DFunc<128>(config);
        } else {
            // raise error
             
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }


    std::chrono::steady_clock::time_point end_time_checkpoint_final = std::chrono::steady_clock::now();
    std::cout << "Time difference (sec) = " <<  
        (std::chrono::duration_cast<std::chrono::microseconds>(end_time_checkpoint_final - begin_time_checkpoint1).count()) /1000000.0  <<std::endl;

    return 0;
}