#pragma once
#include <iostream>
#include <string>
#include <armadillo>
#include <iomanip>
#include <functional>
#include <cmath>

using namespace std;
using namespace arma;

// Analytical 3D Overlap Integral of Primitive Gaussians

// Class that represents a single primitive Gaussian
class Shell {
    public:
        vec center;
        double alpha;
        int angular_momentum;

        Shell(vec center, double alpha, int l);
};

// Class that represents the overlap between two primitive Gaussians
class ShellPair {
    public:
        vec center_a;
        vec center_b;
        vec center_product;
        double alpha;
        double beta;
        int lA;
        int lB;

        // possible angular momentum combinations
        mat lmn_a;
        mat lmn_b;

        // number of functions in each shell
        int func_a;
        int func_b;

        ShellPair(Shell shell_a, Shell shell_b);
        double overlapIntegral1D(int dim, int lA, int lB);
        mat overlapIntegral3D();
};


