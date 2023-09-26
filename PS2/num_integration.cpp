#include "num_integration.h"

/**
 * @brief Constructor for GaussianOverlap class.
 * 
 * This constructor initializes a GaussianOverlap object with the inputted
 * centers, alpha exponents, and angular momentum of two Gaussians.
 * 
 * @param center_a The center of the first Gaussian
 * @param alpha The alpha exponent of the first Gaussian
 * @param lA The angular momentum of the first Gaussian
 * @param center_b The center of the second Gaussian
 * @param beta The alpha exponent of the second Gaussian
 * @param lB The angular momentum of the second Gaussian
 */
GaussianOverlap::GaussianOverlap(double center_a, double alpha, int lA, double center_b, double beta, int lB)
    : center_a(center_a), alpha(alpha), lA(lA), center_b(center_b), beta(beta), lB(lB) {}

/**
 * @brief Calculate 1D Gaussian overlap integrand.
 * 
 * This function calculates the overlap of two Gaussians given their centers, 
 * alpha exponents, and angular momentum.
 * 
 * @param x The coordinate to evaluate the Gaussian overlap at
 * 
 * @return The value of the overlap at x
 */
double GaussianOverlap::overlap1D(double x) {
    double result = pow(x - center_a, lA) * pow(x - center_b, lB) * exp(-alpha * pow(x - center_a, 2) 
                    - beta * pow(x - center_b, 2));

    return result;
}

/**
 * @brief 1D Extended Trapezoidal Rule.
 * 
 * This function performs numerical integration using the extended trapezoidal
 * rule on a given function f(x) over a given interval. 
 * 
 * @param f A function pointer to the function to integrate
 * @param lower_bound The lower bound of the interval
 * @param upper_bound The upper bound of the interval
 * @param num_steps The number of trapezoids to use
 * 
 * @return The value of the integral
 */
double trapezoidal(function<double(double)> f, double lower_bound, double upper_bound, double num_steps) {
    double h = (upper_bound - lower_bound) / num_steps;

    // Factor out the first and last terms
    double integral = (f(lower_bound) + f(upper_bound)) / 2.0;

    // Sum over the rest of the terms
    for (int i = 1; i < num_steps; i++) {
        integral += f(lower_bound + i * h);
    }

    integral *= (upper_bound - lower_bound) / num_steps;
    return integral;
}