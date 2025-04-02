//derivatives.cpp
#include "derivatives.h"


// author of this file: Surya Datta Sudhakar
// department: Computational and Data Sciences
// organization: Indian Institute of Science
// year: 2025
// All the schemes in this header assume uniform grid spacing

// 2 point stencil based
// First order accurate (in error)
// first derivate
// using forward difference scheme
double fds1_first_derivative(double dx, double f_i0, double f_ir1){
    return ((-1 * f_i0) + (1 * f_ir1))/dx;
}

// 2 point stencil based
// First order accurate (in error)
// first derivate
// using backward difference scheme
double bds1_first_derivative(double dx, double f_il1, double f_i0){
    return ((-1 * f_il1) + (1 * f_i0))/dx;
}

// 3 point stencil based
// First order accurate (in error)
// first derivate
// using central difference scheme
double cds2_first_derivative(double dx, double f_il1, double f_i0, double f_ir1){
    return ((-0.5 * f_il1) + (0.5 * f_ir1))/dx;
}

// 3 point stencil based
// Second order accurate (in error)
// first derivate
// using forward difference scheme
double fds2_first_derivative(double dx, double f_i0, double f_ir1, double f_ir2){
    return ((-1.5 * f_i0) + (2 * f_ir1) + (-0.5 * f_ir2))/dx;
}

// 3 point stencil based
// Second order accurate (in error)
// first derivate
// using backward difference scheme
double bds2_first_derivative(double dx, double f_il2, double f_il1, double f_i0){
    return ((0.5 * f_il2) + (-2 * f_il1) + (0.5 * f_i0))/dx;
}

/********************************************************************************/

// 3 point stencil based
// First order accurate (in error)
// second derivate
// using forward difference scheme
double fds1_second_derivative(double dx, double f_i0, double f_ir1, double f_ir2){
    return ((1 * f_i0) + (-2 * f_ir1) + (1 * f_ir2))/(std::pow(dx, 2));
}

// 3 point stencil based
// First order accurate (in error)
// second derivate
// using backward difference scheme
double bds1_second_derivative(double dx, double f_il2, double f_il1, double f_i0){
    return ((1 * f_il2) + (-2 * f_il1) + (1 * f_i0))/(std::pow(dx, 2));
}

// 3 point stencil based
// Second order accurate (in error)
// second derivate
// using central difference scheme
double cds2_second_derivative(double dx, double f_il1, double f_i0, double f_ir1){
    return ((1 * f_il1) + (-2 * f_i0) + (1 * f_ir1))/(std::pow(dx, 2));
}

// 4 point stencil based
// Second order accurate (in error)
// second derivate
// using forward difference scheme
double fds2_second_derivative(double dx, double f_i0, double f_ir1, double f_ir2, double f_ir3){
    return ((2 * f_i0) + (-5 * f_ir1) + (4 * f_ir2) + (-1 * f_ir3))/(std::pow(dx, 2));
}

// 4 point stencil based
// Second order accurate (in error)
// second derivate
// using backward difference scheme
double bds2_second_derivative(double dx, double f_il3, double f_il2, double f_il1, double f_i0){
    return ((-1 * f_il3) + (4 * f_il2) + (-5 * f_il1) + (2 * f_i0))/(std::pow(dx, 2));
}



