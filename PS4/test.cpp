#include <iostream>
#include "molecule.h"
#include "CNDO.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Incorrect file input! << endl";
        return 1;
    }

    const char* filename = argv[1];
    Molecule molecule(filename);
    cout << "Molecule Info: " << endl;
    molecule.printMoleculeInfo();

    cout << "----------------------------------------" << endl;

    // Test all matrix functions before SCF cycle
    cout << "Initial Matrix Calculations: " << endl << endl;
    mat S = calcOverlapMatrix(molecule);
    cout << "Overlap Matrix: " << endl;
    cout << S << endl;
    
    CNDO cndo(molecule, S);
    cout << "Gamma Matrix: " << endl;
    cout << cndo.gammaMatrix << endl;

    cout << "H Core Matrix: " << endl;
    cout << cndo.hCoreMat << endl;

    cout << "Alpha Fock Matrix: " << endl;
    cout << cndo.alphaFockMat << endl;

    cout << "Beta Fock Matrix: " << endl;
    cout << cndo.betaFockMat << endl;

    // Test SCF cycle
    cndo.scfCycle();

    return 0;
}