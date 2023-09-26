#pragma once
#include <iostream>
#include <string>
#include <iomanip>
#include <functional>

using namespace std;

// Numerical Integration: Extended Trapezoidal Rule

// Class that represents the overlap between two Gaussians
class GaussianOverlap {
    public:
        double center_a;
        double center_b;
        double alpha;
        double beta;
        int lA;
        int lB;

        GaussianOverlap(double center_a, double alpha, int lA, double center_b, double beta, int lB);
        double overlap1D(double x);
};

// 1D Extended Trapezoidal Rule
double trapezoidal(function<double(double)> f, double lower_bound, double upper_bound, double num_steps);