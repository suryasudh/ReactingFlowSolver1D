#include "utils_solver.h"

const std::string gas_composition(std::string gas_name){
    if (gas_name == "air"){
        const std::string gas_comp = "H2:0, H:0, O:0, O2:0.20946, OH:0, H2O:0, HO2:0, H2O2:0, AR:0.0097, N2:0.78084";
        return gas_comp;
    } else if (gas_name == "fuel"){
        const std::string gas_comp = "H2:0.41892, H:0, O:0, O2:0.20946, OH:0, H2O:0, HO2:0, H2O2:0, AR:0.0097, N2:0.78084";
        return gas_comp;
    } else if (gas_name == "fuel2"){
        const std::string gas_comp = "H2:2.0, H:0, O:0, O2:1.0, OH:0, H2O:0, HO2:0, H2O2:0, AR:0.0, N2:3.76";
        return gas_comp;
    } else {
        const std::string gas_comp = "H2O:1.0, H2:8.0, AR:1.0";
        return gas_comp;
    }
}



/*
Creates the required initial conditions
// u(x, 0) = sin(2 pi (x/L)) * 0.01
*/
template <typename T>
std::vector<T> initial_cond_on_u_case1(const std::vector<T>& x, T domain_length){
    int grid_size = x.size();
    
    std::vector<T> u(grid_size);

    for (int i = 0; i < grid_size; i++){
        u[i] = static_cast<T>(0.01) * std::sin(2 * 3.14159265358979323846 * (x[i]/domain_length));
    }
    
    return u;
}

/*
Creates the required initial conditions
on temperature with Gaussian distribution 
with peak temperature 1040 K and 400 K at the ends
*/
template <typename T>
std::vector<T> initial_cond_on_temp_case2(const std::vector<T>& x, T domain_length){
    int grid_size = x.size();
    T midpoint_domain = domain_length / static_cast<T>(2);
    T std_dev = static_cast<T>(3.0);

    std::vector<T> temp(grid_size);

    T T_min = static_cast<T>(400);
    T T_peak = static_cast<T>(1200);

    for (int i = 0; i < grid_size; i++){
        T diff = x[i] - midpoint_domain;
        T exponent = - (diff * diff) / (static_cast<T>(2) * std_dev * std_dev);
        T gaussian_value = std::exp(exponent); // Normalized Gaussian (no prefactor needed)

        temp[i] = T_min + (T_peak - T_min) * gaussian_value;
    }

    return temp;
}
