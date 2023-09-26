#include "analytical_integration.h"
#include "factorial.h"

/**
 * @brief Constructor for Shell class.
 * 
 * This constructor initializes a Shell object with the given center,
 * alpha exponent, and angular momentum of the shell.
 * 
 * @param center The center of the shell
 * @param alpha The alpha exponent of the shell
 * @param l The angular momentum of the shell
 */
Shell::Shell(vec center, double alpha, int l) : center(center), alpha(alpha), angular_momentum(l) {}

/**
 * @brief Constructor for ShellPair class.
 * 
 * This constructor initializes a ShellPair object with the information
 * of the two shells. 
 * 
 * @param center_a The center of the first shell
 * @param alpha The alpha exponent of the first shell
 * @param lA The angular momentum of the first shell
 * @param center_b The center of the second shell
 * @param beta The alpha exponent of the second shell
 * @param lB The angular momentum of the second shell
 */
ShellPair::ShellPair(Shell shell_a, Shell shell_b) : center_a(shell_a.center), alpha(shell_a.alpha), lA(shell_a.angular_momentum), 
        center_b(shell_b.center), beta(shell_b.alpha), lB(shell_b.angular_momentum) {

    // Calculate the product of the two centers
    center_product = (alpha * center_a + beta * center_b) / (alpha + beta);

    // Get number of functions in each shell
    if (lA == 0) {
        func_a = 1;
    } else if (lA == 1) {
        func_a = 3;
    } 
    if (lB == 0) {
        func_b = 1;
    } else if (lB == 1) {
        func_b = 3;
    }

    // Get possible angular momentum combinations based on l
    if (lA == 0) {          // s   
        lmn_a = mat({{0, 0, 0}});
    } else if (lA == 1) {   // p 
        lmn_a = mat({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    }
    if (lB == 0) {
        lmn_b = mat({{0, 0, 0}});
    } else if (lB == 1) {
        lmn_b = mat({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    }
}

/**
 * @brief Calculate the analytical integration for the overlap of two primitive Gaussian shells for one dimension.
 * 
 * This function calculates the analytical integration for the overlap of two primitive Gaussian shells
 * for one dimension given the centers, alpha exponents, and angular momentum of the two shells.
 * 
 * @param dim The dimension to calculate the overlap for
 * @param lA The angular momentum of the first shell
 * @param lB The angular momentum of the second shell
 * 
 * @return The value of the overlap
 */
double ShellPair::overlapIntegral1D(int dim, int lA, int lB) {
    // Calculate the exponential prefactor and the associated square root 
    double prefactor = exp(-alpha * beta * pow(center_a(dim) - center_b(dim), 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));

    // Double summation over the angular momentum combinations
    double sum = 0.0;
    for (int i = 0; i <= lA; i++) {
        for (int j = 0; j <= lB; j++) {
            // Only (i + j) even terms contribute
            if ((i + j) % 2 == 0) {
                sum += binomialCoef(lA, i) * binomialCoef(lB, j) * (doubleFactorial(i + j - 1) 
                        * pow(center_product(dim) - center_a(dim), lA - i) * pow(center_product(dim) - center_b(dim), lB - j)) 
                        / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}

/**
 * @brief Calculate the 3D analytical integration for the overlap of two primitive Gaussian shells.
 * 
 * This function evaluates the overlap integral of two primitive Gaussian shells over all dimensions and 
 * over all angular momentum combinations.
 *
 * @return The matrix of overlap integrals
 */
mat ShellPair::overlapIntegral3D() {
    // Initialize matrix of overlap integrals
    mat overlap_integral = zeros(func_a, func_b);

    // Loop over the l, m, n angular momentum combinations of both shells
    for (int i = 0; i < func_a; i++) {
        for (int j = 0; j < func_b; j++) {
            // Calculate overlap integral for each dimension with the current angular momentum combination
            overlap_integral(i, j) = overlapIntegral1D(0, lmn_a(i, 0), lmn_b(j, 0)) * overlapIntegral1D(1, lmn_a(i, 1), lmn_b(j, 1)) 
                    * overlapIntegral1D(2, lmn_a(i, 2), lmn_b(j, 2));
        }
    }
    return overlap_integral;
}