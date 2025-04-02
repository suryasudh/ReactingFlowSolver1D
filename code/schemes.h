//schemes.h
#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include <iomanip>

#ifndef SCHEMES_H
#define SCHEMES_H

// author of this file: Surya Datta Sudhakar
// department: Computational and Data Sciences
// organization: Indian Institute of Science
// year: 2025
// All the schemes in this header assume uniform grid spacing


double explicit_euler(std::function<double (double, double)> func1, double dt, double tn0, double x_tn0);

double implicit_euler(std::function<double (double, double)> func1, double y_tn0, double x_tn0, double dx);

double explicit_rk2(std::function<double (double, double)> func1, double dt, double tn0, double x_tn0);

double explicit_rk4(std::function<double (double, double)> func1, double dt, double tk, double xk);

std::array<double, 2> explicit_rk4_vector(std::function<std::array<double, 2> 
                                                (double, std::array<double, 2>, 
                                                        std::array<double, 5>)> func1, 
                                            double dt, double tn0, 
                                            std::array<double, 2> x_tn0,
                                            std::array<double, 5> props);


double explicit_AB2(std::function<double (double, double)> func1, 
                    double dt, double tn0, double x_tn0
                             , double tnm1, double x_tnm1);

#endif // SCHEMES_H