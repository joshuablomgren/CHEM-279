#pragma once
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include "molecule.h"

using namespace std;
using namespace arma;

// Constants
const double constant_K = 1.75;

// Functions used for symmetric orthogonalization to calculate the energy of the molecule
mat calcHamiltonianMatrix(Molecule molecule, mat overlap_mat);
mat calcXMatrix(mat overlap_mat);
mat calcHamiltonianPrimeMatrix(mat X_mat, mat hamiltonian_mat);
double calcEnergy(mat X_mat, mat hamiltonian_prime_mat, int nElectrons);
double calcEnergy(Molecule molecule, mat overlap_mat);
