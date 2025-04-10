#include <iostream>
#include <vector>       // Although using std::array, vector might be useful for dynamic scenarios if adapted later.
#include <array>        // Primary container for array data.
#include <cmath>        // For std::pow, std::max, M_PI, std::sin
#include <string>       // For std::string
#include <stdexcept>    // For exceptions like std::invalid_argument
#include <numeric>      // For std::accumulate
#include <limits>       // For std::numeric_limits
#include <tuple>        // For returning multiple values from functions

// // --- Placeholder for Cantera Gas Solution Object ---
// // NOTE: This is a very basic placeholder struct. A C++ equivalent of Cantera
// // or a similar thermochemistry/kinetics library (like Cantera's C++ API)
// // would be required for a functional simulation. This struct only mimics
// // the interface used in the Python code.
// template <typename T, size_t N_SPECIES>
// struct GasSolution {
//     // Placeholder properties - these would be calculated based on state in a real library
//     T T_ = 300.0; // Temperature (K)
//     T P_ = 101325.0; // Pressure (Pa)
//     T viscosity_ = 1.8e-5; // Dynamic viscosity (Pa·s)
//     T thermal_conductivity_ = 0.026; // Thermal conductivity (W/(m·K))
//     T mean_molecular_weight_ = 28.97e-3; // Mean molecular weight (kg/mol) - Air approx.
//     std::array<T, N_SPECIES> molecular_weights_; // Molecular weights of species (kg/mol)
//     std::array<T, N_SPECIES> standard_enthalpies_RT_; // Standard molar enthalpies / (R*T) (dimensionless)
//     std::array<T, N_SPECIES> mix_diff_coeffs_; // Mixture-averaged diffusion coefficients (m^2/s) - Approximation
//     std::array<T, N_SPECIES> net_production_rates_; // Net production rates (kmol/(m^3·s))
//     std::array<T, N_SPECIES> X; // Mole fractions (dimensionless) - Should sum to 1

//     // Placeholder for setting state (e.g., based on internal energy U, specific volume V, mass fractions Y)
//     // In Cantera, this would trigger complex calculations.
//     void set_UVY(T internal_energy_mass, T specific_volume, const std::array<T, N_SPECIES>& Y_k) {
//         // --- THIS IS HIGHLY SIMPLIFIED ---
//         // A real implementation would solve for T, P from U, V, Y_k using an equation of state
//         // and then update all thermodynamic and transport properties.
//         // Here, we just assign dummy values for demonstration.
//         T_ = 300.0 + internal_energy_mass / 1000.0; // Fictional update based on energy
//         P_ = GasSolution<T,N_SPECIES>::gas_constant * T_ / (specific_volume * mean_molecular_weight_); // Ideal gas approx.

//         // Update properties based on the new T, P, Y_k (using dummy updates here)
//         viscosity_ = 1.8e-5 * std::pow(T_ / 300.0, 0.7);
//         thermal_conductivity_ = 0.026 * std::pow(T_ / 300.0, 0.8);

//         // Update mole fractions X from mass fractions Y (placeholder calculation)
//         T mass_sum = 0.0;
//         for(size_t k=0; k<N_SPECIES; ++k) {
//              X[k] = Y_k[k] / molecular_weights_[k]; // Proportional to moles
//              mass_sum += X[k]; // Accumulate proportional moles
//         }
//          if (mass_sum > std::numeric_limits<T>::epsilon()) {
//             for(size_t k=0; k<N_SPECIES; ++k) {
//                  X[k] /= mass_sum; // Normalize to get mole fractions
//             }
//          } else {
//              // Handle zero mass case (e.g., uniform distribution)
//              for(size_t k=0; k<N_SPECIES; ++k) X[k] = T{1.0} / N_SPECIES;
//          }


//         // Update other properties like diffusion coeffs, production rates (dummy)
//         mix_diff_coeffs_.fill(1e-5 * std::pow(T_ / 300.0, 1.5));
//         net_production_rates_.fill(0.0); // Assume no reactions for placeholder
//         standard_enthalpies_RT_.fill(0.0); // Assume constant enthalpy base for placeholder
//     }

//     // Accessor methods mirroring Cantera's Python API usage
//     T T() const { return T_; }
//     T P() const { return P_; }
//     T viscosity() const { return viscosity_; }
//     T thermal_conductivity() const { return thermal_conductivity_; }
//     T mean_molecular_weight() const { return mean_molecular_weight_; } // Needs update in set_UVY
//     const std::array<T, N_SPECIES>& molecular_weights() const { return molecular_weights_; }
//     // This method in Cantera returns H_k / (R*T)
//     const std::array<T, N_SPECIES>& standard_enthalpies_RT() const { return standard_enthalpies_RT_; }
//     // This method in Cantera returns mixture-averaged diffusion coeffs D_km
//     const std::array<T, N_SPECIES>& mix_diff_coeffs() const { return mix_diff_coeffs_; }
//     // This method in Cantera returns omega_dot_k in kmol/m^3/s
//     const std::array<T, N_SPECIES>& net_production_rates() const { return net_production_rates_; }

//     // Constructor to initialize arrays (example for air-like mixture)
//     GasSolution() {
//         // Initialize placeholder values (e.g., for N2, O2, Ar...)
//         if constexpr (N_SPECIES > 0) molecular_weights_[0] = 28.0134e-3; // N2
//         if constexpr (N_SPECIES > 1) molecular_weights_[1] = 31.9988e-3; // O2
//         // Fill remaining with placeholder average
//         for(size_t k=2; k<N_SPECIES; ++k) molecular_weights_[k] = 28.97e-3;
//         // Ensure all are initialized if N_SPECIES is small
//         for(size_t k=0; k<N_SPECIES; ++k) {
//             if (molecular_weights_[k] <= 0) molecular_weights_[k] = 28.97e-3;
//         }

//         standard_enthalpies_RT_.fill(0.0);
//         mix_diff_coeffs_.fill(1e-5);
//         net_production_rates_.fill(0.0);
//         // Initial mole fractions (e.g., air-like)
//         if constexpr (N_SPECIES > 0) X[0] = 0.78; // N2
//         if constexpr (N_SPECIES > 1) X[1] = 0.21; // O2
//         T remaining_fraction = 0.01 / (N_SPECIES > 2 ? N_SPECIES - 2 : 1);
//         for(size_t k=2; k<N_SPECIES; ++k) X[k] = remaining_fraction;
//         if constexpr (N_SPECIES == 1) X[0] = 1.0; // Handle single species case

//     }

//     // Universal gas constant J/(mol·K)
//     static constexpr T gas_constant = static_cast<T>(8.31446261815324);
// };


// // --- Helper Functions ---

// // Template function for power, equivalent to np.float_power
// template<typename T>
// T float_power(T base, T exp) {
//     return std::pow(base, exp);
// }

// // Template function for element-wise maximum of an array and a scalar, equivalent to np.maximum(arr, value)
// template<typename T, size_t N>
// std::array<T, N> array_maximum(const std::array<T, N>& arr, T value) {
//     std::array<T, N> result;
//     for (size_t i = 0; i < N; ++i) {
//         result[i] = std::max(arr[i], value);
//     }
//     return result;
// }

// // Template function to sum elements of a std::array, equivalent to np.sum(arr)
// template<typename T, size_t N>
// T array_sum(const std::array<T, N>& arr) {
//     // Start summation with 0.0 of type T
//     return std::accumulate(arr.begin(), arr.end(), static_cast<T>(0.0));
// }


// // --- Converted Core Functions ---

// /**
//  * @brief Calculates the first derivative using central differencing.
//  * @tparam T Floating point type (float, double, long double).
//  * @tparam N Size of the array (number of nodes).
//  * @param arr Input array of values at nodes.
//  * @param dx Grid spacing.
//  * @param bc_type Boundary condition type. Currently only "periodic" is implemented.
//  * @param bc_value Boundary value (not used for periodic).
//  * @return std::array<T, N> containing the calculated derivative at each node.
//  * @throws std::invalid_argument if bc_type is not "periodic".
//  */
// template <typename T, size_t N>
// std::array<T, N> central_diff_first_derivative(
//     const std::array<T, N>& arr,
//     T dx,
//     const std::string& bc_type = "periodic",
//     T bc_value = T{}) // Default constructor for T (0 for numeric types)
// {
//     std::array<T, N> deriv;
//     // Ensure deriv is initialized (equivalent to np.zeros_like)
//     deriv.fill(T{});

//     if constexpr (N == 0) { // Handle empty array case at compile time
//          return deriv;
//     }

//     if (bc_type == "periodic") {
//         // Loop through each node
//         for (size_t i = 0; i < N; ++i) {
//             // Calculate index to the right with wrap-around
//             size_t right = (i + 1) % N;
//             // Calculate index to the left with wrap-around
//             // Equivalent to Python: left = (i - 1) % n
//             size_t left = (i == 0) ? (N - 1) : (i - 1);
//             // Central difference formula
//             deriv[i] = (arr[right] - arr[left]) / (T{2.0} * dx);
//         }
//     }
//     // Add other BC types ('neumann', 'dirichlet') if needed, similar to Python version
//     else {
//         // Throw an error for unsupported boundary condition types
//         throw std::invalid_argument("Invalid bc_type. Only 'periodic' is implemented in this C++ version.");
//     }

//     return deriv;
// }

// /**
//  * @brief Calculates the first derivative using upwind differencing based on velocity sign.
//  * @tparam T Floating point type (float, double, long double).
//  * @tparam N Size of the array (number of nodes).
//  * @param arr Input array of values at nodes.
//  * @param dx Grid spacing.
//  * @param velocity Array of velocities at each node (determines upwind direction).
//  * @param bc_type Boundary condition type. Currently only "periodic" is implemented.
//  * @param bc_value Boundary value (not used for periodic).
//  * @return std::array<T, N> containing the calculated upwind derivative at each node.
//  * @throws std::invalid_argument if bc_type is not "periodic".
//  */
// template <typename T, size_t N>
// std::array<T, N> upwind_first_derivative(
//     const std::array<T, N>& arr,
//     T dx,
//     const std::array<T, N>& velocity,
//     const std::string& bc_type = "periodic",
//     T bc_value = T{})
// {
//     std::array<T, N> arr_ans;
//     // Ensure arr_ans is initialized (equivalent to np.zeros_like)
//     arr_ans.fill(T{});

//     if constexpr (N == 0) { // Handle empty array case at compile time
//          return arr_ans;
//     }

//     if (bc_type == "periodic") {
//         // Loop through each node
//         for (size_t i = 0; i < N; ++i) {
//             // Check velocity direction to determine upwind node
//             if (velocity[i] >= T{0.0}) {
//                 // Positive velocity: Upwind is to the left (i-1)
//                 size_t left = (i == 0) ? (N - 1) : (i - 1); // Periodic BC for left index
//                 arr_ans[i] = (arr[i] - arr[left]) / dx;
//             } else {
//                 // Negative velocity: Upwind is to the right (i+1)
//                 size_t right = (i + 1) % N; // Periodic BC for right index
//                 arr_ans[i] = (arr[right] - arr[i]) / dx;
//             }
//         }
//     } else {
//         // Throw an error for unsupported boundary condition types
//         throw std::invalid_argument("Invalid bc_type. Only 'periodic' is implemented in this C++ version.");
//     }

//     return arr_ans;
// }


// /**
//  * @brief Calculates the time derivatives (dX/dt) for the 1D compressible reacting flow conservation equations.
//  * @tparam T Floating point type (float, double, long double).
//  * @tparam N_NODES Number of spatial grid nodes.
//  * @tparam N_SPECIES Number of chemical species.
//  * @param dx1 Grid spacing.
//  * @param state_cons_old Tuple containing the conserved variables at the current time step:
//  * (rho, rho*u, rho*E, rho*Y_k) at each node.
//  * @param gas1 GasSolution object (placeholder) used to get thermodynamic and transport properties.
//  * Passed by non-const reference as its internal state might be updated (e.g., via set_UVY).
//  * @return std::tuple containing:
//  * - Time derivatives: (d(rho)/dt, d(rho*u)/dt, d(rho*E)/dt, d(rho*Y_k)/dt) at each node.
//  * - Other calculated variables: Tuple(temperature, pressure, production_rates) at each node.
//  */
// template <typename T, size_t N_NODES, size_t N_SPECIES>
// auto conservation_equations_dt(
//     T dx1,
//     const std::tuple<
//         std::array<T, N_NODES>, // rho_old
//         std::array<T, N_NODES>, // rho_u_old
//         std::array<T, N_NODES>, // rho_e0_old (rho * Total Energy E)
//         std::array<std::array<T, N_SPECIES>, N_NODES> // rho_ys_old (rho * Mass Fractions Y_k)
//     >& state_cons_old,
//     GasSolution<T, N_SPECIES>& gas1 // Pass by non-const reference to allow state modification
// ) -> std::tuple< // Return type definition
//         // Time derivatives of conserved variables
//         std::array<T, N_NODES>, // d_rho_dt
//         std::array<T, N_NODES>, // d_rho_u_dt
//         std::array<T, N_NODES>, // d_rho_e0_dt
//         std::array<std::array<T, N_SPECIES>, N_NODES>, // d_rho_ys_dt
//         // Other useful variables calculated during the process
//          std::tuple<
//              std::array<T, N_NODES>, // temperature_old
//              std::array<T, N_NODES>, // pressure_old
//              std::array<std::array<T, N_SPECIES>, N_NODES> // production_rates_old (omega_dot_k * MW_k)
//          > // other_vars tuple
//     >
// {
//     // --- Unpack the input state tuple for easier access ---
//     const auto& [rho_old, rho_u_old, rho_e0_old, rho_ys_old] = state_cons_old;

//     // --- Declare Intermediate Variables (Arrays to store values at each node) ---
//     std::array<T, N_NODES> u_old;               // Velocity u
//     std::array<T, N_NODES> e_int_old;           // Specific internal energy e
//     std::array<T, N_NODES> sp_vol_old;          // Specific volume v = 1/rho

//     std::array<T, N_NODES> temperature_old;     // Temperature T
//     std::array<T, N_NODES> pressure_old;        // Pressure P

//     std::array<T, N_NODES> prop00_mu_old;       // Viscosity mu
//     std::array<T, N_NODES> prop05_th_k_old;     // Thermal conductivity k

//     std::array<T, N_NODES> mean_mix_weight_old; // Mean molecular weight MW_mix
//     std::array<std::array<T, N_SPECIES>, N_NODES> mole_fractions_old; // Mole fractions X_k
//     std::array<std::array<T, N_SPECIES>, N_NODES> mw_species_old;     // Species molecular weights MW_k
//     std::array<std::array<T, N_SPECIES>, N_NODES> enthalpy_k_old;     // Specific enthalpy h_k (J/kg)

//     std::array<std::array<T, N_SPECIES>, N_NODES> mix_diff_coeffs_old;// Mixture diffusion coefficients D_km
//     std::array<std::array<T, N_SPECIES>, N_NODES> production_rates_old;// Production rates omega_dot_k (kmol/m^3/s)


//     // --- Calculate Primitive Variables and Properties at Each Node ---
//     // This loop iterates through each grid point (node)
//     for (size_t i = 0; i < N_NODES; ++i) {
//         // Avoid division by zero if density is extremely small
//         T rho_inv = (rho_old[i] > std::numeric_limits<T>::epsilon()) ? T{1.0} / rho_old[i] : T{0.0};

//         u_old[i] = rho_u_old[i] * rho_inv;                  // u = (rho*u) / rho
//         // Total Energy E = e + u^2/2  =>  e = E - u^2/2
//         e_int_old[i] = (rho_e0_old[i] * rho_inv) - (float_power(u_old[i], T{2.0}) / T{2.0});
//         sp_vol_old[i] = rho_inv;                            // v = 1 / rho

//         // Calculate mass fractions Y_k for the current node
//         std::array<T, N_SPECIES> Y_k_node;
//         for(size_t k=0; k < N_SPECIES; ++k) {
//             // Y_k = (rho*Y_k) / rho
//             Y_k_node[k] = rho_ys_old[i][k] * rho_inv;
//             // Ensure mass fractions are physically plausible (non-negative)
//             Y_k_node[k] = std::max(Y_k_node[k], T{0.0});
//         }
//         // Normalize Y_k_node if needed (due to numerical errors)
//         T Y_sum = array_sum(Y_k_node);
//         if (Y_sum > std::numeric_limits<T>::epsilon()) {
//             for(size_t k=0; k < N_SPECIES; ++k) Y_k_node[k] /= Y_sum;
//         }


//         // Update gas state using the placeholder GasSolution object
//         // This call updates the internal state of gas1 and recalculates properties (T, P, mu, k, etc.)
//         gas1.set_UVY(e_int_old[i], sp_vol_old[i], Y_k_node);

//         // Retrieve properties from the gas object for the current node
//         temperature_old[i] = gas1.T();
//         pressure_old[i] = gas1.P();
//         prop00_mu_old[i] = gas1.viscosity();
//         prop05_th_k_old[i] = gas1.thermal_conductivity();
//         mean_mix_weight_old[i] = gas1.mean_molecular_weight();
//         mole_fractions_old[i] = gas1.X; // Assuming set_UVY updates X based on Y
//         mw_species_old[i] = gas1.molecular_weights(); // Assuming constant or updated by set_UVY
//         mix_diff_coeffs_old[i] = gas1.mix_diff_coeffs(); // Assuming updated by set_UVY
//         production_rates_old[i] = gas1.net_production_rates(); // Assuming updated by set_UVY

//         // Calculate specific enthalpy h_k (J/kg) for each species
//         const auto& h_RT = gas1.standard_enthalpies_RT(); // Get H_k / (R*T)
//         for(size_t k=0; k < N_SPECIES; ++k) {
//              // H_k (J/mol) = (H_k / (R*T)) * R * T
//              T H_k_molar = h_RT[k] * GasSolution<T,N_SPECIES>::gas_constant * temperature_old[i];
//              // h_k (J/kg) = H_k (J/mol) / MW_k (kg/mol)
//              if (mw_species_old[i][k] > std::numeric_limits<T>::epsilon()) { // Avoid division by zero
//                 enthalpy_k_old[i][k] = H_k_molar / mw_species_old[i][k];
//             } else {
//                 enthalpy_k_old[i][k] = T{0.0};
//             }
//         }
//     } // End of loop calculating properties at each node

//     // --- Calculate Spatial Derivatives of Fluxes ---

//     // --- 1. Mass Conservation: d(rho)/dt = - d(rho*u)/dx ---
//     // Calculate gradient of mass flux (rho*u) using central difference
//     auto mass_flux_grad = central_diff_first_derivative(rho_u_old, dx1, "periodic");
//     std::array<T, N_NODES> d_rho_dt;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         d_rho_dt[i] = -mass_flux_grad[i];
//     }

//     // --- 2. Momentum Conservation: d(rho*u)/dt = - d(rho*u*u + P)/dx + d(tau_xx)/dx ---
//     // Convective flux term: rho*u*u
//     std::array<T, N_NODES> convective_flux_mom;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         convective_flux_mom[i] = rho_u_old[i] * u_old[i]; // rho * u * u
//     }
//     // Gradient of convective flux (using upwind based on local velocity u)
//     auto convective_flux_grad_mom = upwind_first_derivative(convective_flux_mom, dx1, u_old, "periodic");

//     // Pressure gradient term: dP/dx (using central difference)
//     auto pressure_grad = central_diff_first_derivative(pressure_old, dx1, "periodic");

//     // Viscous stress term: tau_xx = (4/3)*mu*(du/dx) (assuming Stokes' hypothesis)
//     auto u_grad = central_diff_first_derivative(u_old, dx1, "periodic"); // du/dx
//     std::array<T, N_NODES> viscous_stress_tau_xx; // tau_xx at each node
//     for (size_t i = 0; i < N_NODES; ++i) {
//         viscous_stress_tau_xx[i] = (T{4.0} / T{3.0}) * prop00_mu_old[i] * u_grad[i];
//     }
//     // Gradient of viscous stress term: d(tau_xx)/dx (using central difference)
//     auto viscous_stress_grad_mom = central_diff_first_derivative(viscous_stress_tau_xx, dx1, "periodic");

//     // Combine terms for momentum equation
//     std::array<T, N_NODES> d_rho_u_dt;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         // d(rho*u)/dt = d(tau_xx)/dx - d(rho*u*u)/dx - dP/dx
//         d_rho_u_dt[i] = viscous_stress_grad_mom[i] - convective_flux_grad_mom[i] - pressure_grad[i];
//     }


//     // --- 3. Energy Conservation: d(rho*E)/dt = - d((rho*E + P)*u)/dx + d(u*tau_xx)/dx - d(q)/dx ---
//     // Convective energy flux: (rho*E + P)*u
//     std::array<T, N_NODES> convective_energy_flux;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         convective_energy_flux[i] = (rho_e0_old[i] + pressure_old[i]) * u_old[i];
//     }
//     // Gradient of convective energy flux (using central difference - matching Python code)
//     // Note: Upwind might be more stable for the convective term, but following Python source.
//     auto convective_energy_flux_grad = central_diff_first_derivative(convective_energy_flux, dx1, "periodic");

//     // Viscous work term: u * tau_xx
//     std::array<T, N_NODES> viscous_work_term;
//     for (size_t i = 0; i < N_NODES; ++i) {
//          viscous_work_term[i] = u_old[i] * viscous_stress_tau_xx[i]; // u * tau_xx
//     }
//     // Gradient of viscous work term: d(u*tau_xx)/dx (using central difference)
//     auto viscous_work_grad = central_diff_first_derivative(viscous_work_term, dx1, "periodic");

//     // --- Heat Flux Calculation (q) ---
//     // Heat flux q = q_conduction + q_diffusion = -k*dT/dx + sum(h_k * J_k)
//     // where J_k is the mass diffusion flux of species k.

//     // Calculate temperature gradient: dT/dx (using central difference)
//     auto temp_grad = central_diff_first_derivative(temperature_old, dx1, "periodic");

//     // Calculate conduction heat flux: q_c = -k * dT/dx
//     std::array<T, N_NODES> heat_flux_conduction; // q_c at each node
//     for (size_t i = 0; i < N_NODES; ++i) {
//         heat_flux_conduction[i] = -prop05_th_k_old[i] * temp_grad[i];
//     }

//     // Calculate species concentration gradients: dX_k/dx (using central difference)
//     std::array<std::array<T, N_SPECIES>, N_NODES> mole_fraction_grads; // dXk/dx for each species at each node
//     for (size_t k = 0; k < N_SPECIES; ++k) {
//         std::array<T, N_NODES> Xk_at_nodes; // Temporary array for X_k values across all nodes
//         for (size_t i = 0; i < N_NODES; ++i) {
//             Xk_at_nodes[i] = mole_fractions_old[i][k];
//         }
//         // Calculate gradient for this species
//         auto grad_Xk = central_diff_first_derivative(Xk_at_nodes, dx1, "periodic");
//         // Store the gradient
//         for (size_t i = 0; i < N_NODES; ++i) {
//             mole_fraction_grads[i][k] = grad_Xk[i];
//         }
//     }

//     // Calculate species mass diffusion fluxes J_k
//     // Using Fick's law with mixture-averaged coefficients (approximation): J_k = - rho * D_km * dYk/dx
//     // Or relating to mole fractions: J_k approx - rho * D_km * (MW_k / MW_mix) * dXk/dx
//     // Python code calculates: diff_flux = rho * D * MW * grad(X) / MW_mix
//     // Let's match the Python calculation structure for J_k equivalent term.
//     std::array<std::array<T, N_SPECIES>, N_NODES> diff_flux_species; // Stores the term calculated in Python as 'diff_flux'
//     for (size_t i = 0; i < N_NODES; ++i) {
//         if (mean_mix_weight_old[i] > std::numeric_limits<T>::epsilon()) {
//             T rho_over_mw_mix = rho_old[i] / mean_mix_weight_old[i];
//             for (size_t k = 0; k < N_SPECIES; ++k) {
//                  diff_flux_species[i][k] = rho_over_mw_mix * mix_diff_coeffs_old[i][k] * mw_species_old[i][k] * mole_fraction_grads[i][k];
//             }
//         } else {
//              diff_flux_species[i].fill(T{0.0});
//         }
//     }
//     // Note: The sign convention for J_k might differ. Standard definition often includes a minus sign.
//     // If J_k = - rho * D * ..., then the term calculated above is -J_k (based on standard Fick's law form).
//     // Let's assume diff_flux_species[i][k] represents the J_k used in the Python energy/species equations.

//     // Calculate diffusion heat flux term: q_d = sum(h_k * J_k)
//     std::array<T, N_NODES> diff_flux_enthalpy_sum; // Stores sum(h_k * J_k) at each node
//     diff_flux_enthalpy_sum.fill(T{0.0});
//     for (size_t i = 0; i < N_NODES; ++i) {
//         for (size_t k = 0; k < N_SPECIES; ++k) {
//             // Python code: diff_flux_add = np.sum(enthalpy_k_old * (diff_flux), axis=1)
//             // Assuming diff_flux_species corresponds to Python's diff_flux
//             diff_flux_enthalpy_sum[i] += enthalpy_k_old[i][k] * diff_flux_species[i][k];
//         }
//     }

//     // Total heat flux q = q_c + q_d
//     // Python code: heat_n_diff_flux = ((-1) * heat_flux) + diff_flux_add
//     // where Python's heat_flux = k * dT/dx (missing minus sign for Fourier's law).
//     // So, Python's heat_n_diff_flux = -k*dT/dx + sum(h_k*J_k) = q_c + q_d = q_total
//     std::array<T, N_NODES> heat_flux_total; // q_total at each node
//     for (size_t i = 0; i < N_NODES; ++i) {
//          heat_flux_total[i] = heat_flux_conduction[i] + diff_flux_enthalpy_sum[i];
//     }

//     // Gradient of total heat flux: d(q)/dx (using central difference)
//     auto heat_flux_grad = central_diff_first_derivative(heat_flux_total, dx1, "periodic");

//     // Combine terms for energy equation
//     std::array<T, N_NODES> d_rho_e0_dt;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         // d(rhoE)/dt = - d((rhoE+P)u)/dx + d(u*tau_xx)/dx - d(q)/dx
//         d_rho_e0_dt[i] = -convective_energy_flux_grad[i] + viscous_work_grad[i] - heat_flux_grad[i];
//         // Matching Python signs: d_rho_e0_dt = viscous_work_grad - energy_flux_grad - heat_flux_grad
//         // This implies Python's energy_flux_grad = d((rhoE+P)u)/dx and heat_flux_grad = d(q)/dx
//         // Let's stick to the Python version's signs:
//         // d_rho_e0_dt[i] = viscous_work_grad[i] - convective_energy_flux_grad[i] - heat_flux_grad[i]; // Re-checked, this matches standard form signs if fluxes defined appropriately.
//     }


//     // --- 4. Species Conservation: d(rho*Y_k)/dt = - d(rho*Y_k*u)/dx - d(J_k)/dx + omega_dot_k * MW_k ---
//     // Convective species flux: rho*Y_k*u = (rho*Y_k) * u
//     std::array<std::array<T, N_SPECIES>, N_NODES> convective_species_flux;
//     for(size_t i=0; i < N_NODES; ++i) {
//         for(size_t k=0; k < N_SPECIES; ++k) {
//             convective_species_flux[i][k] = rho_ys_old[i][k] * u_old[i];
//         }
//     }

//     // Gradient of convective species flux: d(rho*Y_k*u)/dx (using central difference - matching Python)
//     // Note: Upwind might be more stable.
//     std::array<std::array<T, N_SPECIES>, N_NODES> convective_species_flux_grad;
//     for (size_t k = 0; k < N_SPECIES; ++k) {
//         std::array<T, N_NODES> flux_k_at_nodes; // Temporary array for flux of species k
//         for (size_t i = 0; i < N_NODES; ++i) {
//             flux_k_at_nodes[i] = convective_species_flux[i][k];
//         }
//         auto grad_flux_k = central_diff_first_derivative(flux_k_at_nodes, dx1, "periodic");
//         for (size_t i = 0; i < N_NODES; ++i) {
//             convective_species_flux_grad[i][k] = grad_flux_k[i];
//         }
//     }

//     // Gradient of diffusion flux: d(J_k)/dx (using central difference)
//     // J_k is represented by diff_flux_species calculated earlier.
//      std::array<std::array<T, N_SPECIES>, N_NODES> diffusion_flux_grad;
//      for (size_t k = 0; k < N_SPECIES; ++k) {
//         std::array<T, N_NODES> Jk_at_nodes; // Diffusion flux J_k (or equivalent term)
//         for (size_t i = 0; i < N_NODES; ++i) {
//             Jk_at_nodes[i] = diff_flux_species[i][k];
//         }
//         auto grad_Jk = central_diff_first_derivative(Jk_at_nodes, dx1, "periodic");
//         for (size_t i = 0; i < N_NODES; ++i) {
//             diffusion_flux_grad[i][k] = grad_Jk[i];
//         }
//     }

//     // Chemical source term: omega_dot_k * MW_k (kmol/m^3/s * kg/mol -> kg/m^3/s)
//     // Note: Cantera returns omega_dot_k in kmol/m^3/s. MW_k is in kg/mol. Need kg/kmol.
//     // MW_k (kg/kmol) = MW_k (kg/mol) * 1000
//     std::array<std::array<T, N_SPECIES>, N_NODES> chemical_source_term; // omega_dot_k * MW_k (in kg/m^3/s)
//      for (size_t i = 0; i < N_NODES; ++i) {
//          for (size_t k = 0; k < N_SPECIES; ++k) {
//             // Python code: prod_rates_mws = production_rates_old * mw_species_old
//             // Assuming production_rates_old is omega_dot_k (kmol/m^3/s) and mw_species_old is MW_k (kg/mol)
//             // Need MW_k in kg/kmol
//             chemical_source_term[i][k] = production_rates_old[i][k] * (mw_species_old[i][k] * T{1000.0});
//          }
//      }

//     // Combine terms for species conservation equation
//     std::array<std::array<T, N_SPECIES>, N_NODES> d_rho_ys_dt;
//     for (size_t i = 0; i < N_NODES; ++i) {
//         for (size_t k = 0; k < N_SPECIES; ++k) {
//              // d(rhoYk)/dt = - d(rhoYk*u)/dx - d(Jk)/dx + omega_dot_k * MW_k
//              // Python: d_rho_ys_dt = -rho_y_u_grad - diff_flux_grad + prod_rates_mws
//              // Assuming Python's rho_y_u_grad = d(rhoYk*u)/dx and diff_flux_grad = d(Jk)/dx
//              d_rho_ys_dt[i][k] = -convective_species_flux_grad[i][k] - diffusion_flux_grad[i][k] + chemical_source_term[i][k];
//         }
//     }

//     // --- Package and Return Results ---
//     // Package the 'other' variables calculated during this step
//     auto other_vars = std::make_tuple(temperature_old, pressure_old, chemical_source_term); // Return source term in kg/m^3/s

//     // Return the time derivatives and the other variables
//     return std::make_tuple(d_rho_dt, d_rho_u_dt, d_rho_e0_dt, d_rho_ys_dt, other_vars);
// }


// /**
//  * @brief Advances the solution state by one time step using a specified numerical scheme.
//  * @tparam T Floating point type (float, double, long double).
//  * @tparam N_NODES Number of spatial grid nodes.
//  * @tparam N_SPECIES Number of chemical species.
//  * @param dt_val Time step size.
//  * @param dx_val Grid spacing.
//  * @param state_cons_old Tuple containing the conserved variables at the beginning of the time step:
//  * (rho, rho*u, rho*E, rho*Y_k) at each node.
//  * @param gas_obj GasSolution object (placeholder) passed to conservation_equations_dt.
//  * @param scheme Time integration scheme ("explicit_euler" or "rk2").
//  * @return std::tuple containing:
//  * - New state: (rho_new, rho_u_new, rho_e0_new, rho_ys_new) at the end of the time step.
//  * - Other variables: Tuple(temperature, pressure, production_rates) corresponding to the *final* derivative evaluation (k1 for Euler, k2 for RK2).
//  * @throws std::invalid_argument if the scheme is not recognized.
//  */
// template <typename T, size_t N_NODES, size_t N_SPECIES>
// auto time_stepper(
//     T dt_val,
//     T dx_val,
//     const std::tuple<
//         std::array<T, N_NODES>, // rho_old
//         std::array<T, N_NODES>, // rho_u_old
//         std::array<T, N_NODES>, // rho_e0_old
//         std::array<std::array<T, N_SPECIES>, N_NODES> // rho_ys_old
//     >& state_cons_old,
//     GasSolution<T, N_SPECIES>& gas_obj, // Pass gas object by reference
//     const std::string& scheme = "rk2"
// ) -> std::tuple< // Return type definition: New state + other vars
//         std::array<T, N_NODES>, // rho_new
//         std::array<T, N_NODES>, // rho_u_new
//         std::array<T, N_NODES>, // rho_e0_new
//         std::array<std::array<T, N_SPECIES>, N_NODES>, // rho_ys_new
//         // Other variables from the final derivative evaluation step
//          std::tuple<
//              std::array<T, N_NODES>, // temperature
//              std::array<T, N_NODES>, // pressure
//              std::array<std::array<T, N_SPECIES>, N_NODES> // production_rates (source term)
//          > // other_vars tuple
//     >
// {
//     // Unpack the initial state
//     const auto& [rho_old, rho_u_old, rho_e0_old, rho_ys_old] = state_cons_old;

//     // --- Calculate k1: Derivatives at the beginning of the step (t = t_n) ---
//     // Call the function that calculates all d(state)/dt terms
//     auto [rho_k1, rho_u_k1, rho_e0_k1, rho_ys_k1, other_vars_k1] =
//         conservation_equations_dt<T, N_NODES, N_SPECIES>(dx_val, state_cons_old, gas_obj);

//     // --- Select Time Integration Scheme ---
//     if (scheme == "explicit_euler") {
//         // --- Explicit Euler Step: y_{n+1} = y_n + dt * f(t_n, y_n) ---
//         std::array<T, N_NODES> rho_new;
//         std::array<T, N_NODES> rho_u_new;
//         std::array<T, N_NODES> rho_e0_new;
//         std::array<std::array<T, N_SPECIES>, N_NODES> rho_ys_new;

//         // Update each conserved variable
//         for(size_t i=0; i<N_NODES; ++i) {
//             rho_new[i] = rho_old[i] + dt_val * rho_k1[i];
//             rho_u_new[i] = rho_u_old[i] + dt_val * rho_u_k1[i];
//             rho_e0_new[i] = rho_e0_old[i] + dt_val * rho_e0_k1[i];
//             for(size_t k=0; k<N_SPECIES; ++k) {
//                 rho_ys_new[i][k] = rho_ys_old[i][k] + dt_val * rho_ys_k1[i][k];
//             }
//         }
//          // Apply density floor (minimum allowed density) - equivalent to np.maximum(rho, 1e-8)
//         rho_new = array_maximum(rho_new, static_cast<T>(1e-8));

//         // Return the new state and the 'other_vars' calculated with k1
//         return std::make_tuple(rho_new, rho_u_new, rho_e0_new, rho_ys_new, other_vars_k1);

//     } else if (scheme == "rk2") {
//         // --- RK2 (Midpoint Method Variant): y* = y_n + dt/2 * k1; k2 = f(t_n+dt/2, y*); y_{n+1} = y_n + dt * k2 ---
//         // This matches the structure of the provided Python code.

//         // --- Calculate Midpoint State (y*) ---
//         std::array<T, N_NODES> rho_mid;
//         std::array<T, N_NODES> rho_u_mid;
//         std::array<T, N_NODES> rho_e0_mid;
//         std::array<std::array<T, N_SPECIES>, N_NODES> rho_ys_mid;

//         T dt_half = dt_val / T{2.0};
//         for(size_t i=0; i<N_NODES; ++i) {
//             rho_mid[i] = rho_old[i] + dt_half * rho_k1[i];
//             rho_u_mid[i] = rho_u_old[i] + dt_half * rho_u_k1[i];
//             rho_e0_mid[i] = rho_e0_old[i] + dt_half * rho_e0_k1[i];
//              for(size_t k=0; k<N_SPECIES; ++k) {
//                 rho_ys_mid[i][k] = rho_ys_old[i][k] + dt_half * rho_ys_k1[i][k];
//             }
//         }
//         // Apply density floor to the midpoint state (as done in Python code)
//         rho_mid = array_maximum(rho_mid, static_cast<T>(1e-8));

//         // Package the midpoint state into a tuple
//         auto state_cons_mid = std::make_tuple(rho_mid, rho_u_mid, rho_e0_mid, rho_ys_mid);

//         // --- Calculate k2: Derivatives at the midpoint state (t = t_n + dt/2) ---
//         // *** IMPORTANT: Using dx_val/2 here, directly following the provided Python code. ***
//         // Standard RK2 would typically use the original dx_val for evaluating k2.
//         // This might be a specific modification or potential error in the original Python.
//         auto [rho_k2, rho_u_k2, rho_e0_k2, rho_ys_k2, other_vars_k2] =
//              conservation_equations_dt<T, N_NODES, N_SPECIES>(dx_val / T{2.0}, state_cons_mid, gas_obj);
//              // Standard RK2 would use: conservation_equations_dt<...>(dx_val, state_cons_mid, gas_obj);

//         // --- Calculate Final State (y_{n+1}) using k2 ---
//         std::array<T, N_NODES> rho_new;
//         std::array<T, N_NODES> rho_u_new;
//         std::array<T, N_NODES> rho_e0_new;
//         std::array<std::array<T, N_SPECIES>, N_NODES> rho_ys_new;

//          for(size_t i=0; i<N_NODES; ++i) {
//             rho_new[i] = rho_old[i] + dt_val * rho_k2[i]; // y_n + dt * k2
//             rho_u_new[i] = rho_u_old[i] + dt_val * rho_u_k2[i];
//             rho_e0_new[i] = rho_e0_old[i] + dt_val * rho_e0_k2[i];
//              for(size_t k=0; k<N_SPECIES; ++k) {
//                 rho_ys_new[i][k] = rho_ys_old[i][k] + dt_val * rho_ys_k2[i][k];
//             }
//         }
//          // Apply density floor to the final state
//         rho_new = array_maximum(rho_new, static_cast<T>(1e-8));

//         // Return the new state and the 'other_vars' calculated with k2
//         return std::make_tuple(rho_new, rho_u_new, rho_e0_new, rho_ys_new, other_vars_k2);

//     } else {
//          // Throw error for unrecognized scheme
//          throw std::invalid_argument("Invalid scheme specified. Choose 'explicit_euler' or 'rk2'.");
//     }
// }



