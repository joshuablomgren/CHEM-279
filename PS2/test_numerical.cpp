#include "num_integration.h"

/**
 * @brief Numerical Integrations 1D overlap integrals of SS and SP Gaussians at origin and offset of 1.
 * 
 * @param a The lower bound of the integration interval
 * @param b The upper bound of the integration interval
 * @param n The number of steps to use
 */
void calcOverlapIntegrals1D(double a, double b, int n) {
    // Create GaussianOverlap objects
    GaussianOverlap ss_origin(0, 1, 0, 0, 1, 0);
    GaussianOverlap sp_origin(0, 1, 0, 0, 1, 1);
    GaussianOverlap ss_offset(0, 1, 0, 1, 1, 0);
    GaussianOverlap sp_offset(0, 1, 0, 1, 1, 1);

    // Bind overlap1D function to pass to trapezoidal function
    auto ss_origin_func = bind(&GaussianOverlap::overlap1D, &ss_origin, placeholders::_1);
    auto sp_origin_func = bind(&GaussianOverlap::overlap1D, &sp_origin, placeholders::_1);
    auto ss_offset_func = bind(&GaussianOverlap::overlap1D, &ss_offset, placeholders::_1);
    auto sp_offset_func = bind(&GaussianOverlap::overlap1D, &sp_offset, placeholders::_1);

    // Calculate overlap integrals using extended trapezoidal rule
    double ss_origin_integral = trapezoidal(ss_origin_func, a, b, n);
    double sp_origin_integral = trapezoidal(sp_origin_func, a, b, n);
    double ss_offset_integral = trapezoidal(ss_offset_func, a, b, n);
    double sp_offset_integral = trapezoidal(sp_offset_func, a, b, n);

    // Print results
    cout << "1D Overlap Integral: a=" << a << ", b=" << b << ", n=" << n << endl;
    cout << "-------------------------------------------------" << endl;
    cout << "- SS Integral at Origin: " << setw(20) << ss_origin_integral << endl;
    cout << "- SP Integral at Origin: " << setw(20) << sp_origin_integral << endl;
    cout << "- SS Integral with Offset of 1: " << setw(13) << ss_offset_integral << endl;
    cout << "- SP Integral with Offset of 1: " << setw(13) << sp_offset_integral << endl;
    cout << endl;
}

int main() {
    cout << "Numerical Integration" << endl;

    // Calculate overlap integrals for different integration limits and number of steps
    calcOverlapIntegrals1D(-5, 5, 5);
    calcOverlapIntegrals1D(-5, 5, 100);
    calcOverlapIntegrals1D(-5, 5, 10000);

    return 0;
}