#include "gold_lj.h"

using namespace std;
using namespace arma;

/**
 * @brief GoldCluster class constructor
 * 
 * This constructor reads in the number of atoms, the atomic numbers of the atoms, and the
 * Cartesian coordinates of the atoms from a file and initializes a GoldCluster object.
 * 
 * @param filename The name of the file containing the cluster data
 */
GoldCluster::GoldCluster(const char* filename) {
    // Open the file
    ifstream infile(filename);

    // Check if the file is valid
    assert(infile.good());

    // Read in the number of atoms
    infile >> num_atoms;

    // Resize the vectors to hold the atomic numbers and coordinates
    z_vals.resize(num_atoms);
    coords.resize(num_atoms, 3);

    // Read in the atomic numbers and coordinates
    for (int i = 0; i < num_atoms; i++) {
        infile >> z_vals[i] >> coords(i, 0) >> coords(i, 1) >> coords(i, 2);
    }

    // Close the file
    infile.close();
}

/**
 * @brief Print the geometry of the cluster
 * 
 * This function prints the number of atoms, the atomic numbers of the atoms, and the
 * Cartesian coordinates of the atoms just as they appear in the input file.
 */
void GoldCluster::print_geometry() {
    cout << "Geometry:" << endl;
    // Print the number of atoms
    cout << num_atoms << endl;

    // Print the atomic numbers and coordinates
    for (int i = 0; i < num_atoms; i++) {
        cout << z_vals[i] << setw(12) << coords(i, 0) << setw(12) << coords(i, 1) 
        << setw(12) << coords(i, 2) << endl;
    }
}

/**
 * @brief Calculate the distance between two atoms
 * 
 * This function calculates the distance between two atoms using their Cartesian coordinates.
 * 
 * @param i The index of the first atom
 * @param j The index of the second atom
 * 
 * @return The distance between the two atoms
 */
double GoldCluster::calcDistance(int i, int j) {
    // Calculate the difference in the x, y, and z coordinates
    double dx = coords(i, 0) - coords(j, 0);
    double dy = coords(i, 1) - coords(j, 1);
    double dz = coords(i, 2) - coords(j, 2);

    // Calculate the distance using the distance formula
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

/**
 * @brief Calculate the Lennard Jones potential between two atoms
 * 
 * This function calculates the Lennard Jones potential between two atoms using their
 * Cartesian coordinates and epsilon and sigma values.
 * 
 * @param i The index of the first atom
 * @param j The index of the second atom
 * 
 * @return The Lennard Jones potential between the two atoms
 */
double GoldCluster::calcLJ(int i, int j) {
    // Calculate the distance between the two atoms
    double dist = calcDistance(i, j);

    // Calculate the Lennard Jones potential
    double lj = EPSILON_AU*(pow(SIGMA_AU/dist, 12) - 2*pow(SIGMA_AU/dist, 6));
    return lj;
}

/**
 * @brief Calculate the total Lennard Jones potential of the cluster
 * 
 * This function calculates the total Lennard Jones potential of the cluster by summing
 * the Lennard Jones potentials between all pairs of atoms.
 * 
 * @return The total Lennard Jones potential of the cluster
 */
double GoldCluster::calcTotalEnergy() {
    // Initialize the total Lennard Jones potential energy
    double energy = 0.0;

    // Loop over all pairs of atoms
    for (int i = 0; i < num_atoms; i++) {
        for (int j = i + 1; j < num_atoms; j++) {
            // Calculate the Lennard Jones potential between the two atoms
            double lj = calcLJ(i, j);

            // Add the Lennard Jones potential to the total
            energy += lj;
        }
    }

    return energy;
}