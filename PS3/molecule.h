#pragma once
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

struct BasisFunction {
    string AO_type;
    rowvec center;
    vec lmn;
    vec exponents;
    vec contraction_coeffs;
    vec norm_constants;

    BasisFunction(string AO_type, rowvec center, vec lmn, vec exponents, vec contraction_coeffs);
    void calcNormConstants();
};

// Class that contains all the information about a molecule
class Molecule {
    public:
        Molecule(const char* filename);
        int nAtoms;
        vec atomicNumbers;
        mat coordinates;   
        int carbon_count;
        int hydrogen_count;
        int nElectrons; 
        int nBasisFunctions;
        vector<BasisFunction> basisFunctionsList;

        // Information on H and C basis functions
        vec H_exponents;
        vec H_coefficients;
        vec C_exponents;
        vec C_2s_coefficients;
        vec C_2p_coefficients;

        void printMoleculeInfo();

    private:
        int countBasisFunctions();
        int countElectrons();
        vector<BasisFunction> buildBasisFunctionsList();
};

// Functions used to calculate overlap integrals
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);
double overlapIntegral3D(rowvec centers_a, rowvec centers_b, double alpha, double beta, vec lmn_a, vec lmn_b);
double calcContractedOverlap(BasisFunction basisFunction1, BasisFunction basisFunction2);
mat calcOverlapMatrix(Molecule molecule);