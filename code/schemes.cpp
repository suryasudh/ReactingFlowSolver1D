#include "schemes.h"
#include "derivatives.h"
#include "utils_solver.h"



/*
Explicit Euler Scheme

--------
Parameters
----------
func1: any function which can be represented as f(t, x_tn0).
dt: the timestep considered. Smaller the timestep, less the error in general.
tn0: state of time in the current step.
x_tn0: state of space in the current step.

Returns
-------
x_tn1

Notes
-----
k1 = func1(tn0, x_tn0)

x_tn1 = x_tn0 + (dt * k1)
*/
double explicit_euler(std::function<double (double, double)> func1, double dt, double tn0, double x_tn0){
    double k1 = func1(tn0, x_tn0);

    double x_tn1 = x_tn0 + (k1 * dt);
    return x_tn1;
}

// Implicit Euler formula: y_{n+1} = y_n + [dx * f(x_{n+1}, y_{n+1})]
// To solve for y_{n+1}, we need to rearrange the equation or use a numerical solver.
double implicit_euler(std::function<double (double, double)> func1, double dx, double x_tn0, double y_tn0) {

    // Initial guess for y_{n+1} (using explicit Euler as a starting point)
    double y_tn1 = y_tn0 + (dx * func1(x_tn0, y_tn0));

    // Tolerance and maximum iterations for the numerical solver
    const double tolerance = 1e-6;
    const int max_iterations = 100;

    // Newton's method to solve for y_{n+1}
    for (int i = 0; i < max_iterations; ++i) {
        // Residual: R = y_tn1 - y_tn0 - dx * f(x_tn0 + dx, y_tn1)
        double residual = y_tn1 - (y_tn0 + (dx * func1(x_tn0 + dx, y_tn1)));

        // Derivative of the residual with respect to y_tn1: dR/dy_tn1 = 1 - dx * df/dy
        // For simplicity, assume df/dy is approximated numerically
        double delta_y = 1e-6; // Small perturbation for numerical derivative
        double df_dy = (func1(x_tn0 + dx, y_tn1 + delta_y) - func1(x_tn0 + dx, y_tn1)) / delta_y;
        double dR_dy = 1 - dx * df_dy;

        // Update y_next using Newton's method
        double y_next_new = y_tn1 - residual / dR_dy;

        // Check for convergence
        if (std::abs(y_next_new - y_tn1) < tolerance) {
            y_tn1 = y_next_new;
            break;
        }

        y_tn1 = y_next_new;
    }

    return y_tn1;
}


/*
2nd Order Runge Kutta Scheme

--------
Parameters
----------
func1: any function which can be represented as f(t, x_tn0).
dt: the timestep considered. Smaller the timestep, less the error in general.
tn0: state of time in the current step.
x_tn0: state of space in the current step.

Returns
-------
x_tn1

Notes
-----
k1 = func1(tn0, x_tn0)
k2 = func1((tn0 + dt), (x_tn0 + (dt * k1)))

x_tn1 = x_tn0 + ((dt/2) * (k1 + k2))
*/
double explicit_rk2(std::function<double (double, double)> func1, double dt, double tn0, double x_tn0){
    double k1 = func1(tn0, x_tn0);
    double k2 = func1((tn0 + dt), (x_tn0 + (k1 * dt)));

    double x_tn1 = x_tn0 + ((k1 + k2) * (dt/2));
    return x_tn1;
}



/*
4th Order Runge Kutta Scheme

--------
Parameters
----------
func1: any function which can be represented as f(t, x_tn0).
dt: the timestep considered. Smaller the timestep, less the error in general.
tn0: state of time in the current step.
x_tn0: state of space in the current step.

Returns
-------
x_tn1

Notes
-----
k1 = func1(tn0, x_tn0)
k2 = func1((tn0 + (dt/2)), (x_tn0 + ((dt/2) * k1)))
k3 = func1((tn0 + (dt/2)), (x_tn0 + ((dt/2) * k2)))
k4 = func1((tn0 + dt), (x_tn0 + (dt * k3)))

x_tn1 = x_tn0 + ((dt/6) * (k1 + (2*k2) + (2*k3) + k4))
*/
double explicit_rk4(std::function<double (double, double)> func1, double dt, double tn0, double x_tn0){
    double k1 = func1(tn0, x_tn0);
    double k2 = func1((tn0 + (dt/2)), (x_tn0 + (k1 * (dt/2))));
    double k3 = func1((tn0 + (dt/2)), (x_tn0 + (k2 * (dt/2))));
    double k4 = func1((tn0 + dt), (x_tn0 + (k3 * dt)));

    double x_tn1 = x_tn0 + ((k1 + (2 * k2) + (2 * k3) + k4) * (dt/6));
    return x_tn1;
}


/*
2nd Order Adam-Bashforth Scheme

--------
Parameters
----------
func1: any function which can be represented as f(t, x_tn0).
dt: the timestep considered. Smaller the timestep, less the error in general.

tn0: state of time in the current step.
x_tn0: state of space in the current step.

tnm1: state of time in the previous step.
x_tnm1: state of space in the previous step.

Returns
-------
x_tn1

Notes
-----
step1 = func1(x_tnm1, tnm1)
step2 = func1(tn0, x_tn0)

x_tn1 = x_tn0 + ((dt/2) * ((3 * step2) - (1 * step1)))
*/
double explicit_AB2(std::function<double (double, double)> func1, 
                    double dt, double tn0, double x_tn0
                             , double tnm1, double x_tnm1){
    double step1 = func1(tnm1, x_tnm1);
    double step2 = func1(tn0, x_tn0);

    double x_tn1 = x_tn0 + ((dt/2) * ((3 * step2) - (1 * step1)));
    return x_tn1;
}




std::array<double, 2> explicit_rk4_vector(std::function<std::array<double, 2> 
                                                (double, std::array<double, 2>, 
                                                         std::array<double, 5>)> func1, 
                                            double dt, 
                                            double tn0, 
                                            std::array<double, 2> x_tn0,
                                            std::array<double, 5> props){

    std::array<double, 2> k1, k2, k3, k4, 
                            xsub1, xsub2, xsub3,
                               x_tn1;

    k1 = func1(tn0, x_tn0, props);
    for (int i = 0; i < xsub1.size(); i++){
        xsub1[i] = x_tn0[i] + (k1[i] * (dt/2));
    }
    
    k2 = func1((tn0 + (dt/2)), (xsub1), props);
    for (int i = 0; i < xsub2.size(); i++){
        xsub2[i] = x_tn0[i] + (k2[i] * (dt/2));
    }

    k3 = func1((tn0 + (dt/2)), (xsub2), props);
    for (int i = 0; i < xsub3.size(); i++){
        xsub3[i] = x_tn0[i] + (k3[i] * dt);
    }

    k4 = func1((tn0 + dt), (xsub3), props);
    for (int i = 0; i < x_tn1.size(); i++){
        x_tn1[i] = x_tn0[i] + ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + k4[i]) * (dt/6));
    }
     
    return x_tn1;
}