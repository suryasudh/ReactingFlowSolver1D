// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include "derivatives.h"
#include "initial_conds.h"
#include "solver_funcs.h"
#include <stdexcept>
#include <type_traits> 

using namespace Cantera;


template <typename T>
class FluidSolver1D {
    static_assert(std::is_floating_point<T>::value, "Template parameter T must be a floating-point type (float, double, long double)");

    public:
        // Member Variables
        std::shared_ptr<Cantera::Solution> gas;
        std::size_t n_species;

        // Time parameters
        T dt;

        // Temperature, Pressure, Mole Fractions of the Gas
        T T0;
        T P0;
        std::string X0;

        // Spatial Discretization
        T L;
        std::size_t Nx;
        T dx;
        std::vector<T> x;
        T max_velocity = static_cast<T>(100.0);
        T density_floor = static_cast<T>(1e-6);

        // Boundary condition mode
        int boundary_mode;

        // Boundary Conditions Dirichlet
        T left_boundary_dirichlet;
        T right_boundary_dirichlet;

        // Boundary Conditions Neumann
        T left_boundary_neumann;
        T right_boundary_neumann;

        std::string scheme_str;
        int current_iteration = 0;

        std::list<std::vector<T>> rho_old_history;
        std::list<std::vector<T>> rho_u_old_history;
        std::list<std::vector<T>> rho_e0_old_history;
        std::list<std::vector<std::vector<T>>> rho_ys_old_history;

        std::tuple<std::list<std::vector<T>>, std::list<std::vector<T>>, 
        std::list<std::vector<T>>, std::list<std::vector<std::vector<T>>>> history_last_five_steps;

        FluidSolver1D(std::shared_ptr<Cantera::Solution>  gas_passed, T T0_passed, T P0_passed, 
                        std::string X0_passed, T L_passed, std::size_t Nx_passed, T dt_passed,
                        std::string scheme_str_passed, int boundary_mode_passed, 
                        T left_boundary_dirichlet_passed, T right_boundary_dirichlet_passed, 
                        T left_boundary_neumann_passed, T right_boundary_neumann_passed) {

            gas = gas_passed;
            n_species = gas->thermo()->nSpecies();
            
            dt = dt_passed;

            T0 = T0_passed;
        
            P0 = P0_passed;
            X0 = X0_passed;
            
            L = L_passed;
            Nx = Nx_passed;
            dx = L / static_cast<T>(Nx);

            x = linspace_utils(static_cast<T>(0.0), L, Nx);  
            
            scheme_str = scheme_str_passed;
            boundary_mode = boundary_mode_passed;
            
            left_boundary_dirichlet = left_boundary_dirichlet_passed;
            right_boundary_dirichlet = right_boundary_dirichlet_passed;
            
            left_boundary_neumann = left_boundary_neumann_passed;
            right_boundary_neumann = right_boundary_neumann_passed;
        }


    

    
        std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
                std::vector<std::vector<T>>> get_initial_vals(int case_num) {

            std::vector<T> q00_rho(Nx);
            std::vector<T> q01_rho_u(Nx);
            std::vector<T> q02_rho_e0(Nx);
            std::vector<std::vector<T>> q03_rho_ys(Nx, std::vector<T>(n_species));
            std::vector<T> sv00_u(Nx);
            std::vector<T> sv03_U(Nx);

            if (case_num == 1){
                const std::string X0_c = X0;

                std::cout << "before setting temp: temp is: " << T0 << " K" << std::endl;
                gas->thermo()->setState_TPX(T0, P0, X0_c);
                std::cout << "after setting temp: temp is: " << gas->thermo()->temperature() << " K" << std::endl;
            
                T initial_density = static_cast<T>(gas->thermo()->density());
                T initial_int_energy = static_cast<T>(gas->thermo()->intEnergy_mass());

                double* ys_double_c = new double[n_species];
                gas->thermo()->getMoleFractions(ys_double_c);
                std::vector<T> initial_ys = vector_typecast<T>(ys_double_c, n_species);

                sv03_U = create_initial_valued_vector(x, initial_int_energy);
                sv00_u = initial_cond_on_u_case1(x, L);

                // std::vector<T> q00_rho(Nx, initial_density);
                q00_rho = create_initial_valued_vector(x, initial_density);
                q01_rho_u = multiply_vectors(q00_rho, sv00_u);
                
                q02_rho_e0 = multiply_vectors(q00_rho, add_vectors(sv03_U,
                    divide_vector_scalar(vector_float_power(sv00_u, static_cast<T>(2.0)), static_cast<T>(2.0))));

                for (int i=0; i<Nx; i++){
                    for (int j=0; j<n_species; j++){
                        q03_rho_ys[i][j] = initial_ys.at(j) * q00_rho.at(i);
                    }
                } 

            } else if (case_num == 2){
                sv00_u = create_initial_valued_vector(q00_rho, static_cast<T>(0.0));
                std::vector<T> temp_init = initial_cond_on_temp_case2(x, L);
                
                // same chemical composition
                // temperature with gaussian peak
                // 0 initial velocity
                // uniform initial pressure
                for (int i=0; i<Nx; i++){
                    const std::string X0_c = X0;
                    gas->thermo()->setState_TPX(temp_init[i], P0, X0_c);
                
                    T initial_density = static_cast<T>(gas->thermo()->density());
                    T initial_int_energy = static_cast<T>(gas->thermo()->intEnergy_mass());

                    double* ys_double_c = new double[n_species];
                    gas->thermo()->getMoleFractions(ys_double_c);
                    std::vector<T> initial_ys = vector_typecast<T>(ys_double_c, n_species);

                    sv03_U[i] = initial_int_energy;
                    

                    // std::vector<T> q00_rho(Nx, initial_density);
                    q00_rho[i] =  initial_density;
                    q01_rho_u[i] = q00_rho[i] * sv00_u[i];

                    q02_rho_e0[i] = q00_rho[i] * (sv03_U[i] + ((sv00_u[i] * sv00_u[i])/static_cast<T>(2.0)));
                    
                    for (int j=0; j<n_species; j++){
                        q03_rho_ys[i][j] = initial_ys.at(j) * q00_rho.at(i);
                    }
                    
                }
            } else if (case_num == 3){

                const std::string X0_c = X0;
                // std::cout << "before setting temp: temp is: " << T0 << " K" << std::endl;
                gas->thermo()->setState_TPX(T0, P0, X0_c);
                // std::cout << "after setting temp: temp is: " << gas->thermo()->temperature() << " K" << std::endl;
            
                T initial_density = static_cast<T>(gas->thermo()->density());
                T initial_int_energy = static_cast<T>(gas->thermo()->intEnergy_mass());

                double* ys_double_c = new double[n_species];
                gas->thermo()->getMoleFractions(ys_double_c);
                std::vector<T> initial_ys = vector_typecast<T>(ys_double_c, n_species);

                sv03_U = create_initial_valued_vector(q00_rho, initial_int_energy);
                sv00_u = create_initial_valued_vector(q00_rho, static_cast<T>(0.0));

                // std::vector<T> q00_rho(Nx, initial_density);
                q00_rho = create_initial_valued_vector(q00_rho, initial_density);
                q01_rho_u = multiply_vectors(q00_rho, sv00_u);
                
                q02_rho_e0 = multiply_vectors(q00_rho, add_vectors(sv03_U,
                    divide_vector_scalar(vector_float_power(sv00_u, static_cast<T>(2)), static_cast<T>(2))));

                for (int i=0; i<Nx; i++){
                    for (int j=0; j<n_species; j++){
                        q03_rho_ys[i][j] = initial_ys.at(j) * q00_rho.at(i);
                    }
                }

            } 
            

            rho_old_history.push_back(q00_rho);
            rho_u_old_history.push_back(q01_rho_u);
            rho_e0_old_history.push_back(q02_rho_e0);
            rho_ys_old_history.push_back(q03_rho_ys);

            return std::make_tuple(q00_rho, q01_rho_u, q02_rho_e0, q03_rho_ys);
        }


    
        std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
            std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> solver_func(const std::vector<T>& q00_rho, 
                const std::vector<T>& q01_rho_u, const std::vector<T>& q02_rho_e0, const std::vector<std::vector<T>>& q03_rho_ys) {
            current_iteration += 1;

            // std::cout << "entering solver_func setting temp: temp is: " << gas->thermo()->temperature() << " K" << std::endl;

            std::vector<T> rho_old = q00_rho;
            std::vector<T> rho_u_old = q01_rho_u;
            std::vector<T> rho_e0_old = q02_rho_e0;
            std::vector<std::vector<T>> rho_ys_old = q03_rho_ys;

            std::vector<T> rho_new;
            std::vector<T> rho_u_new;
            std::vector<T> rho_e0_new;
            std::vector<std::vector<T>> rho_ys_new;
            std::vector<T> temperature_old;
            std::vector<T> pressure_old;

            history_last_five_steps = std::make_tuple(rho_old_history, rho_u_old_history, rho_e0_old_history, rho_ys_old_history);
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old, pressure_old) = time_stepper(dt, dx, rho_old, rho_u_old, rho_e0_old, rho_ys_old, gas, scheme_str, current_iteration, history_last_five_steps);

            // std::cout << "after coming back from time_stepper in solver_func setting temp: temp is: " << gas->thermo()->temperature() << " K" << std::endl;
            
            if (current_iteration < 5){
                rho_old_history.push_back(rho_new);
                rho_u_old_history.push_back(rho_u_new);
                rho_e0_old_history.push_back(rho_e0_new);
                rho_ys_old_history.push_back(rho_ys_new);
            } else {
                rho_old_history.pop_front();
                rho_u_old_history.pop_front();
                rho_e0_old_history.pop_front();
                rho_ys_old_history.pop_front();

                rho_old_history.push_back(rho_new);
                rho_u_old_history.push_back(rho_u_new);
                rho_e0_old_history.push_back(rho_e0_new);
                rho_ys_old_history.push_back(rho_ys_new);
            }
            
            return std::make_tuple(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old, pressure_old);
        }


    private:
    // Optional: Add private helper functions if needed for calculations within the class

};

