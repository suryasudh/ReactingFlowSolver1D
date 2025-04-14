#include "utils_solver.h"

const std::string gas_composition(std::string gas_name){
    if (gas_name == "air"){
        const std::string gas_comp = "H2:0, H:0, O:0, O2:0.20946, OH:0, H2O:0, HO2:0, H2O2:0, AR:0.0097, N2:0.78084";
        return gas_comp;
    } else if (gas_name == "fuel"){
        const std::string gas_comp = "H2:0.41892, H:0, O:0, O2:0.20946, OH:0, H2O:0, HO2:0, H2O2:0, AR:0.0097, N2:0.78084";
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
std::vector<T> initial_cond_on_u_case0(std::vector<T> x, T domain_length){
    int grid_size = x.size() - 1;
    
    std::vector<double> u(grid_size+1);

    for (int i = 0; i < grid_size+1; i++){
        u[i] = static_cast<T>(0);
    }
    
    return u;
}


/*
Creates the required initial conditions
// u(x, 0) = sin(2 pi (x/L)) * 0.01
*/
template <typename T>
std::vector<T> initial_cond_on_u_case1(std::vector<T> x, T domain_length){
    int grid_size = x.size() - 1;
    
    std::vector<double> u(grid_size+1);

    for (int i = 0; i < grid_size+1; i++){
        u[i] = static_cast<T>(0.01) * std::sin(2 * 3.14159265358979323846 * (x[i]/domain_length));
    }
    
    return u;
}
