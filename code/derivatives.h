// derivatives.h
#pragma

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include <iomanip>
#ifndef DERIVATIVES_H
#define DERIVATIVES_H

// author of this file: Surya Datta Sudhakar
// department: Computational and Data Sciences
// organization: Indian Institute of Science
// year: 2025

// ========================================================================== //
//                       FIRST DERIVATIVE APPROXIMATIONS                      //
// ========================================================================== //

/**
 * @brief Calculates the first derivative using a 2-point forward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_i0 Function value at the current point f(x).
 * @param f_ir1 Function value at the next point f(x+h).
 * @return First derivative approximation (Order O(h)).
 */
template <typename T>
T fds1_first_derivative(T dx, T f_i0, T f_ir1) {
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    return (f_ir1 - f_i0) / dx; // Simplified (-1*f_i0 + 1*f_ir1)
}

/**
 * @brief Calculates the first derivative using a 2-point backward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x).
 * @return First derivative approximation (Order O(h)).
 */
template <typename T>
T bds1_first_derivative(T dx, T f_il1, T f_i0) {
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    return (f_i0 - f_il1) / dx; // Simplified (-1*f_il1 + 1*f_i0)
}

/**
 * @brief Calculates the first derivative using a 3-point central difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x). (Note: Not used in standard formula)
 * @param f_ir1 Function value at the next point f(x+h).
 * @return First derivative approximation (Order O(h^2)).
 */
template <typename T>
T cds2_first_derivative(T dx, T f_il1, T f_i0, T f_ir1) {
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    // f_i0 is unused in the standard central difference formula for the first derivative.
    return (f_ir1 - f_il1) / (T(2.0) * dx); // Simplified (0.5*f_ir1 - 0.5*f_il1) / dx
}

/**
 * @brief Calculates the first derivative using a 3-point forward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_i0 Function value at the current point f(x).
 * @param f_ir1 Function value at the next point f(x+h).
 * @param f_ir2 Function value at the point f(x+2h).
 * @return First derivative approximation (Order O(h^2)).
 */
template <typename T>
T fds2_first_derivative(T dx, T f_i0, T f_ir1, T f_ir2) {
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    return (T(-1.5) * f_i0 + T(2.0) * f_ir1 - T(0.5) * f_ir2) / dx;
    // Note: Can also be written as (-3*f_i0 + 4*f_ir1 - f_ir2) / (2*dx)
}

/**
 * @brief Calculates the first derivative using a 3-point backward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il2 Function value at the point f(x-2h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x).
 * @return First derivative approximation (Order O(h^2)).
 */
template <typename T>
T bds2_first_derivative(T dx, T f_il2, T f_il1, T f_i0) {
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    return (T(0.5) * f_il2 - T(2.0) * f_il1 + T(1.5) * f_i0) / dx;
}


// ========================================================================== //
//                      SECOND DERIVATIVE APPROXIMATIONS                      //
// ========================================================================== //

/**
 * @brief Calculates the second derivative using a 3-point forward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_i0 Function value at the current point f(x).
 * @param f_ir1 Function value at the next point f(x+h).
 * @param f_ir2 Function value at the point f(x+2h).
 * @return Second derivative approximation (Order O(h)).
 */
template <typename T>
T fds1_second_derivative(T dx, T f_i0, T f_ir1, T f_ir2){
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    T dx_squared = dx * dx;
    return (f_i0 - T(2.0) * f_ir1 + f_ir2) / dx_squared; // Simplified
}

/**
 * @brief Calculates the second derivative using a 3-point backward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il2 Function value at the point f(x-2h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x).
 * @return Second derivative approximation (Order O(h)).
 */
template <typename T>
T bds1_second_derivative(T dx, T f_il2, T f_il1, T f_i0){
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    T dx_squared = dx * dx;
    return (f_il2 - T(2.0) * f_il1 + f_i0) / dx_squared; // Simplified
}

/**
 * @brief Calculates the second derivative using a 3-point central difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x).
 * @param f_ir1 Function value at the next point f(x+h).
 * @return Second derivative approximation (Order O(h^2)).
 */
template <typename T>
T cds2_second_derivative(T dx, T f_il1, T f_i0, T f_ir1){
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    T dx_squared = dx * dx;
    return (f_il1 - T(2.0) * f_i0 + f_ir1) / dx_squared; // Simplified
}

/**
 * @brief Calculates the second derivative using a 4-point forward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_i0 Function value at the current point f(x).
 * @param f_ir1 Function value at the next point f(x+h).
 * @param f_ir2 Function value at the point f(x+2h).
 * @param f_ir3 Function value at the point f(x+3h).
 * @return Second derivative approximation (Order O(h^2)).
 */
template <typename T>
T fds2_second_derivative(T dx, T f_i0, T f_ir1, T f_ir2, T f_ir3){
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    T dx_squared = dx * dx;
    return (T(2.0) * f_i0 - T(5.0) * f_ir1 + T(4.0) * f_ir2 - T(1.0) * f_ir3) / dx_squared;
}

/**
 * @brief Calculates the second derivative using a 4-point backward difference scheme.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param dx Step size (h).
 * @param f_il3 Function value at the point f(x-3h).
 * @param f_il2 Function value at the point f(x-2h).
 * @param f_il1 Function value at the previous point f(x-h).
 * @param f_i0 Function value at the current point f(x).
 * @return Second derivative approximation (Order O(h^2)).
 */
template <typename T>
T bds2_second_derivative(T dx, T f_il3, T f_il2, T f_il1, T f_i0){
    static_assert(std::is_floating_point<T>::value || std::is_integral<T>::value,
                  "Template parameter T must be a numeric type.");
    T dx_squared = dx * dx;
    return (T(-1.0) * f_il3 + T(4.0) * f_il2 - T(5.0) * f_il1 + T(2.0) * f_i0) / dx_squared;
}



// ========================================================================== //
//              COMPUTING 1D DOMAIN DERIVATIVES APPROXIMATIONS                //
// ========================================================================== //



/**
 * @brief Computes first derivatives using 2nd-order central differences with periodic BCs
 * @tparam T Floating-point type (float, double, long double)
 * @param arr Input array of function values
 * @param n Number of elements in the array
 * @param dx Grid spacing
 * @return T* Newly allocated array containing derivatives (caller must manage memory)
 */
template <typename T>
T* central_diff_first_derivative(const T* arr, int n, T dx, int mode) {
    static_assert(std::is_floating_point<T>::value,
                 "Template parameter T must be a floating-point type");

    if(n <= 0) throw std::invalid_argument("Array size must be positive");
    if(dx == T(0)) throw std::invalid_argument("dx cannot be zero");

    T* deriv = new T[n];  // Allocate new array

    if (mode == 1) {
        for(int i = 0; i < n; ++i) {
            const int right = (i + 1) % n;
            const int left = (i - 1 + n) % n;  // Handle negative indices
            deriv[i] = (arr[right] - arr[left]) / (T(2) * dx);
        }
    }

    return deriv;
}




#endif // DERIVATIVES_H