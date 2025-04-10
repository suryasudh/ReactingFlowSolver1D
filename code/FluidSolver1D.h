// main.cpp
#include "cantera/core.h"
#include "utils_solver.h"
#include "derivatives.h"
#include "initial_conds.h"
#include "schemes.h"
#include <stdexcept>
#include <type_traits> 

using namespace Cantera;


template <typename T>
class FluidSolver1D {
    static_assert(std::is_floating_point<T>::value, "Template parameter T must be a floating-point type (float, double, long double)");

    public:
        // Member Variables
        std::shared_ptr<Cantera::ThermoPhase> gas; 
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

        std::string scheme_str;

        FluidSolver1D(std::shared_ptr<Cantera::ThermoPhase> gas_passed, T T0_passed, T P0_passed, 
                        std::string X0_passed, T L_passed, std::size_t Nx_passed, T dt_passed,
                        std::string scheme_str_passed) {
            gas = gas_passed;
            n_species = gas->nSpecies();
            
            dt = dt_passed;

            T0 = T0_passed;
            P0 = P0_passed;
            X0 = X0_passed;
            
            L = L_passed;
            Nx = Nx_passed;
            dx = L / static_cast<T>(Nx);

            x = linspace_utils(static_cast<T>(0.0), L, Nx);    
            
            scheme_str = scheme_str_passed;
        }


    

    
        std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, std::vector<std::vector<T>>> get_initial_vals(int case_num) {
            const std::string X0_c = X0;
            gas->setState_TPX(T0, P0, X0_c);
        

            T initial_density = static_cast<T>(gas->density());
            T initial_int_energy = static_cast<T>(gas->intEnergy_mass());

            

            double* ys_double_c = new double[n_species];
            gas->getMoleFractions(ys_double_c);
            std::vector<T> initial_ys = vector_typecast<T>(ys_double_c, n_species);
            delete[] ys_double_c; 

            std::vector<T> sv03_U(Nx, initial_int_energy);
            std::vector<T> sv00_u = initial_cond_on_u_case1(x, L);

            std::vector<T> q00_rho(Nx, initial_density);
            std::vector<T> q01_rho_u = multiply_vectors(q00_rho, sv00_u);
            
            std::vector<T> q02_rho_e0 = multiply_vectors(q00_rho, add_vectors(sv03_U,
                divide_vector_scalar(vector_float_power(sv00_u, static_cast<T>(2)), static_cast<T>(2))));

            
            
            std::vector<std::vector<T>> q03_rho_ys(Nx, std::vector<T>(n_species));
            for (int i=0; i<Nx; i++){
                for (int j=0; j<n_species; j++){
                    q03_rho_ys[i][j] = initial_ys.at(j) * q00_rho.at(i);
                }
            }


            return std::make_tuple(q00_rho, q01_rho_u, q02_rho_e0, q03_rho_ys);
        }


    
        std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
                    std::vector<std::vector<T>>, std::vector<T>> solver_func(std::vector<T> q00_rho, std::vector<T> q01_rho_u, 
                                                                            std::vector<T> q02_rho_e0, std::vector<std::vector<T>> q03_rho_ys) {
            std::vector<T> rho_old = q00_rho;
            std::vector<T> rho_u_old = q01_rho_u;
            std::vector<T> rho_e0_old = q02_rho_e0;
            std::vector<std::vector<T>> rho_ys_old = q03_rho_ys;

            std::vector<T> rho_new;
            std::vector<T> rho_u_new;
            std::vector<T> rho_e0_new;
            std::vector<std::vector<T>> rho_ys_new;
            std::vector<T> temperature_old;

            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old) = time_stepper(dt, dx, 
                                                                                                    rho_old, rho_u_old,
                                                                                                    rho_e0_old, rho_ys_old, gas, scheme_str);
            
            return std::make_tuple(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old);
        }


    private:
    // Optional: Add private helper functions if needed for calculations within the class

};

