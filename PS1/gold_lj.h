#pragma once
#include <vector>
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>

using namespace std;
using namespace arma;          

#define EPSILON_AU 5.29 * 4.184          // Epsilon of gold converted from kcal/mol to kJ/mol
#define SIGMA_AU 2.951                   // Sigma of gold in angstroms


// Create a class for the cluster of gold atoms
class GoldCluster {
    public:
        int num_atoms;                   // Number of atoms in the cluster
        vector<int> z_vals;              // Atomic numbers of the atoms
        mat coords;                      // Cartesian coordinates of the atoms

        // Constructor
        GoldCluster(const char* filename);

        // Member functions
        void print_geometry();
        double calcDistance(int i, int j);
        double calcLJ(int i, int j);
        double calcTotalEnergy();
};