#include <iostream>
#include <vector>       
#include <deque>
#include <array>        // Primary container for array data.
#include <cmath>        // For std::pow, std::max, M_PI, std::sin
#include <string>       // For std::string
#include <stdexcept>    // For exceptions like std::invalid_argument
#include <numeric>      // For std::accumulate
#include <limits>       // For std::numeric_limits
#include <tuple>        // For returning multiple values from functions
#include "utils_solver.h"
#include "derivatives.h"
#include <cantera/core.h>
#include "cantera/base/AnyMap.h"
#include "cantera/transport/Transport.h"



template <typename T>
std::tuple<int, int, std::vector<T>, std::vector<T>, std::vector<T>, 
std::vector<T>, std::vector<T>, std::vector<T>, std::vector<T>,
std::vector<T>, std::vector<std::vector<T>>, std::vector<std::vector<T>>, std::vector<std::vector<T>>,
std::vector<std::vector<T>>, std::vector<std::vector<T>>> setStateGetProps(const std::vector<T>& rho_old, 
                                                            const std::vector<T>& rho_u_old, 
                                                            const std::vector<T>& rho_e0_old, 
                                                            const std::vector<std::vector<T>>& rho_ys_old, 
                                                            std::shared_ptr<Cantera::Solution> gas_obj1){
    
    
    int n_nodes = rho_old.size();
    int n_species1 = gas_obj1->thermo()->nSpecies();

    std::vector<T> u_old = divide_vectors(rho_u_old, rho_old);
    std::vector<T> e_int_old = subtract_vectors(divide_vectors(rho_e0_old, rho_old),
                        divide_vector_scalar(vector_float_power(u_old, static_cast<T>(2.0)), static_cast<T>(2.0)));
    std::vector<T> sp_vol_old = divide_scalar_vector(static_cast<T>(1.0), rho_old);
    
    std::vector<T> temperature_old(n_nodes, static_cast<T>(0.0));
    std::vector<T> pressure_old(n_nodes, static_cast<T>(0.0));
    
    std::vector<T> prop00_mu_old(n_nodes, static_cast<T>(0.0));
    std::vector<T> prop05_th_k_old(n_nodes, static_cast<T>(0.0));
    
    std::vector<T> mean_mix_weight_old(n_nodes, static_cast<T>(0.0));
    std::vector<std::vector<T>> mole_fractions_old(n_nodes, std::vector<T>(n_species1));
    std::vector<std::vector<T>> mw_species_old(n_nodes, std::vector<T>(n_species1));
    std::vector<std::vector<T>> enthalpy_k_old(n_nodes, std::vector<T>(n_species1));
    
    std::vector<std::vector<T>> mix_diff_coeffs_old(n_nodes, std::vector<T>(n_species1));
    std::vector<std::vector<T>> production_rates_old(n_nodes, std::vector<T>(n_species1));
    
    
    for (int i = 0; i < n_nodes; ++i) {
        std::vector<T> ys_kth_node_old = divide_vector_scalar(rho_ys_old[i], rho_old[i]);
        int _n_size;
        double* _ys_kth_node_old_c;
        std::tie(_ys_kth_node_old_c, _n_size) = vector_typecast_inv(ys_kth_node_old);

        double u_value = static_cast<double>(e_int_old[i]);
        double v_value = static_cast<double>(sp_vol_old[i]);
        gas_obj1->thermo()->setMassFractions(_ys_kth_node_old_c);
        gas_obj1->thermo()->setState_UV(u_value, v_value, 1.0E-9);
        
        
        temperature_old[i] = static_cast<T>(gas_obj1->thermo()->temperature());
        pressure_old[i] = static_cast<T>(gas_obj1->thermo()->pressure());
        
        prop00_mu_old[i] = gas_obj1->transport()->viscosity();
        prop05_th_k_old[i] = gas_obj1->transport()->thermalConductivity();

        mean_mix_weight_old[i] = gas_obj1->thermo()->meanMolecularWeight();

        double* yMoleFractions_double_c = new double[n_species1];
        gas_obj1->thermo()->getMoleFractions(yMoleFractions_double_c);
        mole_fractions_old[i] = vector_typecast<T>(yMoleFractions_double_c, n_species1);

        std::vector<double> yMolecularWeights_double_c = gas_obj1->thermo()->molecularWeights();
        mw_species_old[i] = vector_vector_typecast<T>(yMolecularWeights_double_c);

        double* yEnthalpyRT_double_c = new double[n_species1];
        gas_obj1->thermo()->getEnthalpy_RT(yEnthalpyRT_double_c);
        enthalpy_k_old[i] = vector_typecast<T>(yEnthalpyRT_double_c, n_species1);
        
        double* yMixDiffCoeffs_double_c = new double[n_species1];
        gas_obj1->transport()->getMixDiffCoeffsMass(yMixDiffCoeffs_double_c);
        mix_diff_coeffs_old[i] = vector_typecast<T>(yMixDiffCoeffs_double_c, n_species1);

        double* yProductionRates_double_c = new double[n_species1];
        gas_obj1->kinetics()->getNetProductionRates(yProductionRates_double_c);
        production_rates_old[i] = vector_typecast<T>(yProductionRates_double_c, n_species1);
    }

    return std::make_tuple(n_nodes, n_species1, u_old, e_int_old, sp_vol_old, 
        temperature_old, pressure_old, prop00_mu_old, prop05_th_k_old,
        mean_mix_weight_old, mole_fractions_old, mw_species_old, enthalpy_k_old,
        mix_diff_coeffs_old, production_rates_old);
}


template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> conservation_equations_dt(T dx1, const std::vector<T>& rho_old, 
        const std::vector<T>& rho_u_old, const std::vector<T>& rho_e0_old, 
        const std::vector<std::vector<T>>& rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj1){

    // Setting the properties
    int n_nodes, n_species1;
    std::vector<T> u_old, e_int_old, sp_vol_old;
    std::vector<T> temperature_old, pressure_old;
    std::vector<T> prop00_mu_old, prop05_th_k_old;
    std::vector<T> mean_mix_weight_old;
    std::vector<std::vector<T>> mole_fractions_old, mw_species_old, enthalpy_k_old;
    std::vector<std::vector<T>> mix_diff_coeffs_old, production_rates_old;

    std::tie(n_nodes, n_species1, u_old, e_int_old, sp_vol_old, 
        temperature_old, pressure_old, prop00_mu_old, prop05_th_k_old,
        mean_mix_weight_old, mole_fractions_old, mw_species_old, enthalpy_k_old,
        mix_diff_coeffs_old, production_rates_old) = setStateGetProps(rho_old, rho_u_old, rho_e0_old, rho_ys_old, gas_obj1);    
    

    // Conservation equations
    std::vector<T> d_rho_dt(n_nodes, static_cast<T>(0.0));
    std::vector<T> d_rho_u_dt(n_nodes, static_cast<T>(0.0));
    std::vector<T> d_rho_e0_dt(n_nodes, static_cast<T>(0.0));
    std::vector<std::vector<T>> d_rho_ys_dt(n_nodes, std::vector<T>(n_species1, static_cast<T>(0.0)));

    /****************************************** 
     *** MASS CONSERVATION
    ******************************************/
    std::vector<T> mass_flux_grad = arr_first_derivative_cds_periodic(rho_u_old, dx1);
    d_rho_dt = multiply_vector_scalar(mass_flux_grad, static_cast<T>(-1.0));




    /****************************************** 
     *** MOMENTUM CONSERVATION
    ******************************************/
    std::vector<T> convective_flux = multiply_vectors(rho_u_old, u_old);
    std::vector<T> convective_flux_grad = arr_upwind_first_derivative_periodic(convective_flux, u_old, dx1);

    std::vector<T> pressure_grad = arr_first_derivative_cds_periodic(pressure_old, dx1);

    std::vector<T> viscous_flux = multiply_vector_scalar(
        multiply_vectors(arr_first_derivative_cds_periodic(u_old, dx1), prop00_mu_old), 
        (static_cast<T>(4.0) / static_cast<T>(3.0)));
    std::vector<T> viscous_flux_grad = arr_first_derivative_cds_periodic(viscous_flux, dx1);

    d_rho_u_dt = subtract_vectors((viscous_flux_grad, convective_flux_grad), pressure_grad);
    
    
    /****************************************** 
     *** ENERGY CONSERVATION
    ******************************************/
    std::vector<T> energy_flux = multiply_vectors(add_vectors(rho_e0_old, pressure_old), u_old);
    std::vector<T> energy_flux_grad = arr_first_derivative_cds_periodic(energy_flux, dx1);


    std::vector<T> viscous_work = multiply_vectors(viscous_flux, u_old);
    
    std::vector<T> viscous_work_grad = arr_first_derivative_cds_periodic(viscous_work, dx1);
    std::vector<std::vector <T>> conc_grad = arr_first_derivative_cds_periodic(mole_fractions_old, dx1);

    
    std::vector<std::vector <T>> diff_flux = divide_vectors_2nd_order(multiply_vectors_2nd_order(multiply_vectors_2nd_order(
                                multiply_vectors_2nd_order(tiling_vector(rho_old, n_species1), mix_diff_coeffs_old),
                                mw_species_old), conc_grad), tiling_vector(mean_mix_weight_old, n_species1));
    std::vector<T> diff_flux_add = vector_reduction_by_sum1(multiply_vectors_2nd_order(enthalpy_k_old, diff_flux));
    

    std::vector<T> heat_flux = multiply_vectors(prop05_th_k_old, arr_first_derivative_cds_periodic(temperature_old, dx1));
    std::vector<T> heat_n_diff_flux = add_vectors(multiply_vector_scalar(heat_flux, static_cast<T>(-1.0)), diff_flux_add);

    std::vector<T> heat_flux_grad = arr_first_derivative_cds_periodic(heat_n_diff_flux, dx1);

    d_rho_e0_dt = subtract_vectors(subtract_vectors(viscous_work_grad, energy_flux_grad), heat_flux_grad);
    


    /****************************************** 
     *** SPECIES CONSERVATION
    ******************************************/       
    std::vector<std::vector<T>> rho_y_u_old = multiply_vectors_2nd_order(rho_ys_old, tiling_vector(u_old, n_species1));
    std::vector<std::vector<T>> rho_y_u_grad = arr_first_derivative_cds_periodic(rho_y_u_old, dx1);

    std::vector<std::vector<T>> diff_flux_grad = arr_first_derivative_cds_periodic(diff_flux, dx1);
    std::vector<std::vector<T>> prod_rates_mws = multiply_vectors_2nd_order(production_rates_old, mw_species_old);

    d_rho_ys_dt = add_vectors(subtract_vectors(
                    multiply_vector_2nd_order_scalar(rho_y_u_grad, static_cast<T>(-1.0)), diff_flux_grad), prod_rates_mws);
    
    
    return std::make_tuple(d_rho_dt, d_rho_u_dt, d_rho_e0_dt, d_rho_ys_dt, temperature_old, pressure_old);
}






// Implementing the Runge-Kutta 1st Order accurate method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeRK1(T dt_val, T dx_val, const std::vector<T>& rho_old, const std::vector<T>& rho_u_old, 
        const std::vector<T>& rho_e0_old, const std::vector<std::vector<T>>& rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj){

    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    std::vector<T> rho_k1;
    std::vector<T> rho_u_k1;
    std::vector<T> rho_e0_k1;
    std::vector<std::vector<T>> rho_ys_k1;
    std::vector<T> temperature_old_1;
    std::vector<T> pressure_old_1;

    std::tie(rho_k1, rho_u_k1, rho_e0_k1, rho_ys_k1, 
        temperature_old_1, pressure_old_1) = conservation_equations_dt(dx_val, rho_old, rho_u_old, rho_e0_old, 
                                                                                rho_ys_old, gas_obj);
    
    _rho_new = add_vectors(rho_old, multiply_vector_scalar(rho_k1, dt_val));
    _rho_u_new = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k1, dt_val));
    _rho_e0_new = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k1, dt_val));
    _rho_ys_new = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k1, dt_val));
    _temperature_old = temperature_old_1;
    _pressure_old = pressure_old_1;
    

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}


// Implementing the Runge-Kutta 2nd Order accurate method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeRK2(T dt_val, T dx_val, const std::vector<T>& rho_old, const std::vector<T>& rho_u_old, 
        const std::vector<T>& rho_e0_old, const std::vector<std::vector<T>>& rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj){

    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    std::vector<T> rho_k1;
    std::vector<T> rho_u_k1;
    std::vector<T> rho_e0_k1;
    std::vector<std::vector<T>> rho_ys_k1;
    std::vector<T> temperature_old_1;
    std::vector<T> pressure_old_1;

    std::tie(rho_k1, rho_u_k1, rho_e0_k1, rho_ys_k1, temperature_old_1, pressure_old_1) = conservation_equations_dt(dx_val, rho_old, rho_u_old, rho_e0_old, rho_ys_old, gas_obj);

    std::vector<T> rho_dt_k1;
    std::vector<T> rho_u_dt_k1;
    std::vector<T> rho_e0_dt_k1;
    std::vector<std::vector<T>> rho_ys_dt_k1;
    std::vector<T> temperature_old_dt_k1;
    std::vector<T> pressure_old_dt_k1;


    // Values at 1/2 step in time
    T half_timestep = dt_val/static_cast<T>(2.0);
    rho_dt_k1 = add_vectors(rho_old, multiply_vector_scalar(rho_k1, half_timestep));
    rho_u_dt_k1 = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k1, half_timestep));
    rho_e0_dt_k1 = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k1, half_timestep));
    rho_ys_dt_k1 = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k1, half_timestep));
    temperature_old_dt_k1 = temperature_old_1;
    pressure_old_dt_k1 = pressure_old_1;

    // Determining the slopes at the half step in time
    std::vector<T> rho_k2;
    std::vector<T> rho_u_k2;
    std::vector<T> rho_e0_k2;
    std::vector<std::vector<T>> rho_ys_k2;
    std::vector<T> temperature_old_2;
    std::vector<T> pressure_old_2;

    std::tie(rho_k2, rho_u_k2, rho_e0_k2, rho_ys_k2,                   // made a change here to the dx_val by removing the /2 part
        temperature_old_2, pressure_old_2) = conservation_equations_dt(dx_val, rho_dt_k1, rho_u_dt_k1, rho_e0_dt_k1, rho_ys_dt_k1, gas_obj);
    
    // Taking a step now from starting position using the slopes at the half step
    _rho_new = add_vectors(rho_old, multiply_vector_scalar(rho_k2, dt_val));
    _rho_u_new = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k2, dt_val));
    _rho_e0_new = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k2, dt_val));
    _rho_ys_new = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k2, dt_val));
    _temperature_old = temperature_old_1;
    _pressure_old = pressure_old_1;
    

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}



// Implementing the Runge-Kutta 3rd Order accurate method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeRK3(T dt_val, T dx_val, const std::vector<T>& rho_old, const std::vector<T>& rho_u_old, 
        const std::vector<T>& rho_e0_old, const std::vector<std::vector<T>>& rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj){

    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    std::vector<T> rho_k1;
    std::vector<T> rho_u_k1;
    std::vector<T> rho_e0_k1;
    std::vector<std::vector<T>> rho_ys_k1;
    std::vector<T> temperature_old_1;
    std::vector<T> pressure_old_1;

    std::tie(rho_k1, rho_u_k1, rho_e0_k1, rho_ys_k1, temperature_old_1, pressure_old_1) = conservation_equations_dt(dx_val, rho_old, rho_u_old, rho_e0_old, rho_ys_old, gas_obj);

    std::vector<T> rho_dt_k1;
    std::vector<T> rho_u_dt_k1;
    std::vector<T> rho_e0_dt_k1;
    std::vector<std::vector<T>> rho_ys_dt_k1;
    std::vector<T> temperature_old_dt_k1;
    std::vector<T> pressure_old_dt_k1;

    // Values at 1/2 step in time
    // y0 + k1/2
    T half_timestep = dt_val/static_cast<T>(2.0);
    rho_dt_k1 = add_vectors(rho_old, multiply_vector_scalar(rho_k1, half_timestep));
    rho_u_dt_k1 = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k1, half_timestep));
    rho_e0_dt_k1 = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k1, half_timestep));
    rho_ys_dt_k1 = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k1, half_timestep));
    temperature_old_dt_k1 = temperature_old_1;
    pressure_old_dt_k1 = pressure_old_1;

    // Determining the slopes at the half step in time
    std::vector<T> rho_k2;
    std::vector<T> rho_u_k2;
    std::vector<T> rho_e0_k2;
    std::vector<std::vector<T>> rho_ys_k2;
    std::vector<T> temperature_old_2;
    std::vector<T> pressure_old_2;

    // k2 without dt i.e. still need to multiply dt_val or its multiple
    std::tie(rho_k2, rho_u_k2, rho_e0_k2, rho_ys_k2, 
        temperature_old_2, pressure_old_2) = conservation_equations_dt(dx_val, rho_dt_k1, rho_u_dt_k1, rho_e0_dt_k1, rho_ys_dt_k1, gas_obj);
    
    std::vector<T> rho_dt_k2;
    std::vector<T> rho_u_dt_k2;
    std::vector<T> rho_e0_dt_k2;
    std::vector<std::vector<T>> rho_ys_dt_k2;
    std::vector<T> temperature_old_dt_k2;
    std::vector<T> pressure_old_dt_k2;

    // Values at 1/2 step in time
    // y0 + (2k2 - k1)
    rho_dt_k2 = add_vectors(rho_old, subtract_vectors(multiply_vector_scalar(rho_k2, dt_val*static_cast<T>(2.0)), multiply_vector_scalar(rho_k1, dt_val)));
    rho_u_dt_k2 = add_vectors(rho_u_old, subtract_vectors(multiply_vector_scalar(rho_u_k2, dt_val*static_cast<T>(2.0)), multiply_vector_scalar(rho_u_k1, dt_val)));
    rho_e0_dt_k2 = add_vectors(rho_e0_old, subtract_vectors(multiply_vector_scalar(rho_e0_k2, dt_val*static_cast<T>(2.0)), multiply_vector_scalar(rho_e0_k1, dt_val)));
    rho_ys_dt_k2 = add_vectors(rho_ys_old, subtract_vectors(multiply_vector_2nd_order_scalar(rho_ys_k2, dt_val*static_cast<T>(2.0)), multiply_vector_2nd_order_scalar(rho_ys_k1, dt_val)));
    temperature_old_dt_k2 = temperature_old_2;
    pressure_old_dt_k2 = pressure_old_2;

    // Determining the slopes at the half step in time
    std::vector<T> rho_k3;
    std::vector<T> rho_u_k3;
    std::vector<T> rho_e0_k3;
    std::vector<std::vector<T>> rho_ys_k3;
    std::vector<T> temperature_old_3;
    std::vector<T> pressure_old_3;

    std::tie(rho_k3, rho_u_k3, rho_e0_k3, rho_ys_k3, 
        temperature_old_3, pressure_old_3) = conservation_equations_dt(dx_val, rho_dt_k2, rho_u_dt_k2, rho_e0_dt_k2, rho_ys_dt_k2, gas_obj);
    

    // Taking a step now from starting position using the slopes at the half step
    _rho_new = add_vectors(rho_old, add_vectors(add_vectors(multiply_vector_scalar(rho_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_k2, dt_val/static_cast<T>(1.5))),
                                                                        multiply_vector_scalar(rho_k3, dt_val/static_cast<T>(6.0))));
    _rho_u_new = add_vectors(rho_u_old, add_vectors(add_vectors(multiply_vector_scalar(rho_u_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_u_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_u_k3, dt_val/static_cast<T>(3.0))));
    _rho_e0_new = add_vectors(rho_e0_old, add_vectors(add_vectors(multiply_vector_scalar(rho_e0_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_e0_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_e0_k3, dt_val/static_cast<T>(3.0))));
    _rho_ys_new = add_vectors(rho_ys_old, add_vectors(add_vectors(multiply_vector_2nd_order_scalar(rho_ys_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_2nd_order_scalar(rho_ys_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_2nd_order_scalar(rho_ys_k3, dt_val/static_cast<T>(3.0))));
    _temperature_old = temperature_old_1;
    _pressure_old = pressure_old_1;  

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}



// Implementing the Runge-Kutta 4th Order accurate method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeRK4(T dt_val, T dx_val, const std::vector<T>& rho_old, const std::vector<T>& rho_u_old, 
        const std::vector<T>& rho_e0_old, const std::vector<std::vector<T>>& rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj){

    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;
    
    std::vector<T> rho_k1;
    std::vector<T> rho_u_k1;
    std::vector<T> rho_e0_k1;
    std::vector<std::vector<T>> rho_ys_k1;
    std::vector<T> temperature_old_1;
    std::vector<T> pressure_old_1;

    std::tie(rho_k1, rho_u_k1, rho_e0_k1, rho_ys_k1, temperature_old_1, pressure_old_1) = conservation_equations_dt(dx_val, rho_old, rho_u_old, rho_e0_old, rho_ys_old, gas_obj);

    std::vector<T> rho_dt_k1;
    std::vector<T> rho_u_dt_k1;
    std::vector<T> rho_e0_dt_k1;
    std::vector<std::vector<T>> rho_ys_dt_k1;
    std::vector<T> temperature_old_dt_k1;
    std::vector<T> pressure_old_dt_k1;

    // Values at 1/2 step in time
    T half_timestep = dt_val/static_cast<T>(2.0);
    rho_dt_k1 = add_vectors(rho_old, multiply_vector_scalar(rho_k1, half_timestep));
    rho_u_dt_k1 = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k1, half_timestep));
    rho_e0_dt_k1 = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k1, half_timestep));
    rho_ys_dt_k1 = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k1, half_timestep));
    temperature_old_dt_k1 = temperature_old_1;
    pressure_old_dt_k1 = pressure_old_1;

    // Determining the slopes at the half step in time
    std::vector<T> rho_k2;
    std::vector<T> rho_u_k2;
    std::vector<T> rho_e0_k2;
    std::vector<std::vector<T>> rho_ys_k2;
    std::vector<T> temperature_old_2;
    std::vector<T> pressure_old_2;

    std::tie(rho_k2, rho_u_k2, rho_e0_k2, rho_ys_k2, 
        temperature_old_2, pressure_old_2) = conservation_equations_dt(dx_val, rho_dt_k1, rho_u_dt_k1, rho_e0_dt_k1, rho_ys_dt_k1, gas_obj);
    
    std::vector<T> rho_dt_k2;
    std::vector<T> rho_u_dt_k2;
    std::vector<T> rho_e0_dt_k2;
    std::vector<std::vector<T>> rho_ys_dt_k2;
    std::vector<T> temperature_old_dt_k2;
    std::vector<T> pressure_old_dt_k2;

    // Values at 1/2 step in time
    rho_dt_k2 = add_vectors(rho_old, multiply_vector_scalar(rho_k2, half_timestep));
    rho_u_dt_k2 = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k2, half_timestep));
    rho_e0_dt_k2 = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k2, half_timestep));
    rho_ys_dt_k2 = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k2, half_timestep));
    temperature_old_dt_k2 = temperature_old_2;
    pressure_old_dt_k2 = pressure_old_2;

    // Determining the slopes at the half step in time
    std::vector<T> rho_k3;
    std::vector<T> rho_u_k3;
    std::vector<T> rho_e0_k3;
    std::vector<std::vector<T>> rho_ys_k3;
    std::vector<T> temperature_old_3;
    std::vector<T> pressure_old_3;

    std::tie(rho_k3, rho_u_k3, rho_e0_k3, rho_ys_k3, 
        temperature_old_3, pressure_old_3) = conservation_equations_dt(dx_val, rho_dt_k2, rho_u_dt_k2, rho_e0_dt_k2, rho_ys_dt_k2, gas_obj);
    
    std::vector<T> rho_dt_k3;
    std::vector<T> rho_u_dt_k3;
    std::vector<T> rho_e0_dt_k3;
    std::vector<std::vector<T>> rho_ys_dt_k3;
    std::vector<T> temperature_old_dt_k3;
    std::vector<T> pressure_old_dt_k3;

    // Values at 1 step in time
    rho_dt_k3 = add_vectors(rho_old, multiply_vector_scalar(rho_k3, dt_val));
    rho_u_dt_k3 = add_vectors(rho_u_old, multiply_vector_scalar(rho_u_k3, dt_val));
    rho_e0_dt_k3 = add_vectors(rho_e0_old, multiply_vector_scalar(rho_e0_k3, dt_val));
    rho_ys_dt_k3 = add_vectors(rho_ys_old, multiply_vector_2nd_order_scalar(rho_ys_k3, dt_val));
    temperature_old_dt_k3 = temperature_old_3;
    pressure_old_dt_k3 = pressure_old_3;
                                                                        
    // Determining the slopes at the 1 step in time
    std::vector<T> rho_k4;
    std::vector<T> rho_u_k4;
    std::vector<T> rho_e0_k4;
    std::vector<std::vector<T>> rho_ys_k4;
    std::vector<T> temperature_old_4;
    std::vector<T> pressure_old_4;

    std::tie(rho_k4, rho_u_k4, rho_e0_k4, rho_ys_k4, 
        temperature_old_4, pressure_old_4) = conservation_equations_dt(dx_val, rho_dt_k3, rho_u_dt_k3, rho_e0_dt_k3, rho_ys_dt_k3, gas_obj); 
    

    // Taking a step now from starting position using the slopes at the half step
    _rho_new = add_vectors(rho_old, add_vectors(add_vectors(add_vectors(multiply_vector_scalar(rho_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_k3, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_k4, dt_val/static_cast<T>(6.0))));
    _rho_u_new = add_vectors(rho_u_old, add_vectors(add_vectors(add_vectors(multiply_vector_scalar(rho_u_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_u_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_u_k3, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_u_k4, dt_val/static_cast<T>(6.0))));
    _rho_e0_new = add_vectors(rho_e0_old, add_vectors(add_vectors(add_vectors(multiply_vector_scalar(rho_e0_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_scalar(rho_e0_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_e0_k3, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_scalar(rho_e0_k4, dt_val/static_cast<T>(6.0))));
    _rho_ys_new = add_vectors(rho_ys_old, add_vectors(add_vectors(add_vectors(multiply_vector_2nd_order_scalar(rho_ys_k1, dt_val/static_cast<T>(6.0)),
                                                                        multiply_vector_2nd_order_scalar(rho_ys_k2, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_2nd_order_scalar(rho_ys_k3, dt_val/static_cast<T>(3.0))),
                                                                        multiply_vector_2nd_order_scalar(rho_ys_k4, dt_val/static_cast<T>(6.0))));
    _temperature_old = temperature_old_4;
    _pressure_old = pressure_old_4;

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}


// Implementing the Adams-Bashforth method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeAB2(T dt_val, T dx_val, std::vector<T> rho_old, std::vector<T> rho_u_old, 
    std::vector<T> rho_e0_old, std::vector<std::vector<T>> rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj,
    std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, 
    std::vector<std::vector<T>>, std::vector<std::vector<std::vector<T>>>> history_1st_five_steps){
    
    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}


// Implementing the Adams-Bashforth method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeAB3(T dt_val, T dx_val, std::vector<T> rho_old, std::vector<T> rho_u_old, 
    std::vector<T> rho_e0_old, std::vector<std::vector<T>> rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj,
    std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, 
    std::vector<std::vector<T>>, std::vector<std::vector<std::vector<T>>>> history_1st_five_steps){
    
    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}


// Implementing the Adams-Bashforth method
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> schemeAB4(T dt_val, T dx_val, std::vector<T> rho_old, std::vector<T> rho_u_old, 
    std::vector<T> rho_e0_old, std::vector<std::vector<T>> rho_ys_old, std::shared_ptr<Cantera::Solution> gas_obj,
    std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, 
    std::vector<std::vector<T>>, std::vector<std::vector<std::vector<T>>>> history_1st_five_steps){
    
    std::vector<T> _rho_new;
    std::vector<T> _rho_u_new;
    std::vector<T> _rho_e0_new;
    std::vector<std::vector<T>> _rho_ys_new;
    std::vector<T> _temperature_old;
    std::vector<T> _pressure_old;

    _rho_new = vector_value_minimum_limiter(_rho_new, static_cast<T>(1e-8));
    _rho_ys_new = vector_value_minimum_limiter(_rho_ys_new, static_cast<T>(0.0));

    return std::make_tuple(_rho_new, _rho_u_new, _rho_e0_new, _rho_ys_new, _temperature_old, _pressure_old);
}






template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, 
    std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> time_stepper(T dt_val, T dx_val, std::vector<T> rho_old, std::vector<T> rho_u_old,
                                                                                std::vector<T> rho_e0_old, std::vector<std::vector<T>> rho_ys_old,
                                                                                std::shared_ptr<Cantera::Solution> gas_obj, std::string scheme, int current_iteration,
                                                                                std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, 
                                                                                std::vector<std::vector<T>>, std::vector<std::vector<std::vector<T>>>> history_1st_five_steps){
    std::vector<T> rho_new;
    std::vector<T> rho_u_new;
    std::vector<T> rho_e0_new;
    std::vector<std::vector<T>> rho_ys_new;
    std::vector<T> temperature_old;
    std::vector<T> pressure_old;

    if (scheme == "RK1"){
        std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
            temperature_old, pressure_old) = schemeRK1(dt_val, dx_val, rho_old, rho_u_old, 
                rho_e0_old, rho_ys_old, gas_obj);

    } else if (scheme == "RK2"){
        std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
            temperature_old, pressure_old) = schemeRK2(dt_val, dx_val, rho_old, rho_u_old, 
                rho_e0_old, rho_ys_old, gas_obj);  

    } else if (scheme == "RK3") {
        std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
            temperature_old, pressure_old) = schemeRK3(dt_val, dx_val, rho_old, rho_u_old, 
                rho_e0_old, rho_ys_old, gas_obj);  

    } else if (scheme == "RK4") {
        std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
            temperature_old, pressure_old) = schemeRK4(dt_val, dx_val, rho_old, rho_u_old, 
                rho_e0_old, rho_ys_old, gas_obj);         

    } else if (scheme == "AB2"){
        if (current_iteration < 2){
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeRK2(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj); 
        } else {
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeAB2(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj, history_1st_five_steps); 
        }

    } else if (scheme == "AB3"){
        if (current_iteration < 3){
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeRK3(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj); 
        } else {
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeAB3(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj, history_1st_five_steps); 
        }

    } else if (scheme == "AB4"){
        if (current_iteration < 4){
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeRK4(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj); 
        } else {
            std::tie(rho_new, rho_u_new, rho_e0_new, rho_ys_new, 
                temperature_old, pressure_old) = schemeAB4(dt_val, dx_val, rho_old, rho_u_old, 
                    rho_e0_old, rho_ys_old, gas_obj, history_1st_five_steps); 
        }

    } else {
        // Throw error for unrecognized scheme
        throw std::invalid_argument("Invalid scheme specified. Choose 'explicit_euler' or 'rk2'.");
    }
    

    return std::make_tuple(rho_new, rho_u_new, rho_e0_new, rho_ys_new, temperature_old, pressure_old);              
}









