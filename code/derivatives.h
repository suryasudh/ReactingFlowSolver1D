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

// 2 point stencil based
// First order accurate (in error)
// first derivate
// using forward difference scheme
double fds1_first_derivative(double dx, double f_i0, double f_ir1);

// 2 point stencil based
// First order accurate (in error)
// first derivate
// using backward difference scheme
double bds1_first_derivative(double dx, double f_il1, double f_i0);

// 3 point stencil based
// First order accurate (in error)
// first derivate
// using central difference scheme
double cds2_first_derivative(double dx, double f_il1, double f_i0, double f_ir1);

// 3 point stencil based
// Second order accurate (in error)
// first derivate
// using forward difference scheme
double fds2_first_derivative(double dx, double f_i0, double f_ir1, double f_ir2);

// 3 point stencil based
// Second order accurate (in error)
// first derivate
// using backward difference scheme
double bds2_first_derivative(double dx, double f_il2, double f_il1, double f_i0);

/********************************************************************************/

// 3 point stencil based
// First order accurate (in error)
// second derivate
// using forward difference scheme
double fds1_second_derivative(double dx, double f_i0, double f_ir1, double f_ir2);

// 3 point stencil based
// First order accurate (in error)
// second derivate
// using backward difference scheme
double bds1_second_derivative(double dx, double f_il2, double f_il1, double f_i0);

// 3 point stencil based
// Second order accurate (in error)
// second derivate
// using central difference scheme
double cds2_second_derivative(double dx, double f_il1, double f_i0, double f_ir1);

// 4 point stencil based
// Second order accurate (in error)
// second derivate
// using forward difference scheme
double fds2_second_derivative(double dx, double f_i0, double f_ir1, double f_ir2, double f_ir3);

// 4 point stencil based
// Second order accurate (in error)
// second derivate
// using backward difference scheme
double bds2_second_derivative(double dx, double f_il3, double f_il2, double f_il1, double f_i0);




#endif // DERIVATIVES_H