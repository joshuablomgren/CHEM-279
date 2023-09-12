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
 * Cartesian coordinates of the atoms.
 */
void GoldCluster::print_geometry() {
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
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
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
    double lj = EPSILON_AU * (pow(SIGMA_AU / dist, 12) - 2 * pow(SIGMA_AU / dist, 6));
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

    // Loop over all the atom pairs
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

/**
 * @brief Calculate the force between two atoms using analytical gradient
 * 
 * This function calculates the force between two atoms using their Cartesian coordinates
 * and epsilon and sigma values and the analytical gradient of the Lennard Jones potential
 * with respect to the x, y, and z coordinates.
 * 
 * @param i The index of the first atom
 * @param j The index of the second atom
 * 
 * @return The force between the two atoms
 */
vec GoldCluster::calcForce(int i, int j) {
    // Calculate the distance between the two atoms
    double dist = calcDistance(i, j);

    // Calculate the force using the gradient of the Lennard Jones potential
    double lj_grad = -12 * EPSILON_AU * (pow(SIGMA_AU, 12) / pow(dist, 13) - pow(SIGMA_AU, 6) / pow(dist, 7));

    // Calculate the derivative of the distance formula with respect to the x, y, and z coordinates
    // Note: This is positive for atom i and negative for atom j, so this will be added for atom i and subtracted for atom j
    // in the force matrix calculation
    double dDistX = (coords(i, 0) - coords(j, 0)) / dist;
    double dDistY = (coords(i, 1) - coords(j, 1)) / dist;
    double dDistZ = (coords(i, 2) - coords(j, 2)) / dist;

    // Calculate the force components
    double force_x = -lj_grad * dDistX;
    double force_y = -lj_grad * dDistY;
    double force_z = -lj_grad * dDistZ;

    // Return the force vector
    vec force_vec = {force_x, force_y, force_z};
    return force_vec;
}

/**
 * @brief Calculate the force matrix of the cluster
 * 
 * This function calculates the force matrix of the cluster by calculating the force
 * between all pairs of atoms i and j. 
 * 
 * @return The force matrix of the cluster
 */
mat GoldCluster::calcAnalyticalForces() {
    // Initialize the force matrix
    mat force_matrix(num_atoms, 3);

    // Loop over all the atom pairs
    for (int i = 0; i < num_atoms; i++) {
        for (int j = i + 1; j < num_atoms; j++) {
            // Calculate the force between the two atoms
            vec force_vec = calcForce(i, j);

            // Force vector is added for atom i and subtracted for atom j because the force on atom i is
            // opposite to the force on atom j
            // Transpose to add force vector to correct coord column
            force_matrix.row(i) += force_vec.t();   
            force_matrix.row(j) -= force_vec.t();
        }
    }
    return force_matrix;
}

/**
 * @brief Calculate the force matrix using forward difference approximation
 * 
 * This function calculates the force matrix of the cluster by calculating the force
 * on each atom using the forward difference approximation.
 *
 * @param step The step size for the forward difference approximation
 * 
 * @return The force matrix of the cluster approximated using forward difference
 */
mat GoldCluster::forwardDifference(double step) {
    mat force_matrix(num_atoms, 3);

    // Get the initial energy E(x)
    double E_initial = calcTotalEnergy();

    // Loop over all the atoms
    for (int i = 0; i < num_atoms; i++) {
        // Loop over all the coordinates
        for (int j = 0; j < 3; j++) {
            // Add the step size to the coordinate
            coords(i, j) += step;

            // Calculate energy with new coordinate E(x + h)
            double E_forward = calcTotalEnergy();

            // Calculate force using the forward difference approximation
            double force = -(E_forward - E_initial) / step;
            force_matrix(i, j) = force;

            // Set coordinate back
            coords(i, j) -= step;
        }
    }
    return force_matrix;
}

/**
 * @brief Calculate the force matrix using central difference approximation
 * 
 * This function calculates the force matrix of the cluster by calculating the force
 * on each atom using the central difference approximation.
 * 
 * @param step The step size for the central difference approximation
 * 
 * @return The force matrix of the cluster approximated using central difference
 */
mat GoldCluster::centralDifference(double step) {
    mat force_matrix(num_atoms, 3);

    double E_initial = calcTotalEnergy();

    for (int i = 0; i < num_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            coords(i, j) += step;

            // E(x + h)
            double E_forward = calcTotalEnergy();

            // Set coordinate back by two steps
            coords(i, j) -= 2 * step;

            // E(x - h)
            double E_backward = calcTotalEnergy();

            // Calculate force using the central difference approximation
            double force = -(E_forward - E_backward) / (2 * step);  // negative gradient
            force_matrix(i, j) = force;

            // Set coordinate back
            coords(i, j) += step;
        }
    }
    return force_matrix;
}

/**
 * @brief Steepest descent with Golden Section Line Search
 * 
 * This function performs steepest descent with Golden Section Line Search to find the
 * minimum energy configuration of the cluster.
 * 
 * @param step The intiial step size for the steepest descent
 * @param tol The tolerance for the steepest descent
 * 
 * @return The minimum energy configuration of the cluster
 */
mat GoldCluster::steepestDescent(double stepsize, double tol) {
    // Count iterations
    int iter = 0;

    // Evaluate the initial energy and force matrix
    double old_energy = calcTotalEnergy();
    double new_energy;
    mat force_matrix = calcAnalyticalForces();

    // Save the old coordinates
    mat old_coords;

    // Print initial geometry, energy, and norm of the force matrix
    cout << "Initial Geometry before Steepest Descent:" << endl;
    print_geometry();
    cout << "Initial LJ Energy: " << old_energy << endl;
    cout << "Initial Norm of Forces: " << norm(force_matrix, "fro") << endl << endl;

    // Steepest descent loop
    while (norm(force_matrix, "fro") > tol and iter < 1000) {
        // Save the old coordinates
        old_coords = coords;
        // Calculate the new coordinates
        mat new_coords = coords + stepsize * force_matrix / norm(force_matrix);

        // Set the coordinates to the new coordinates
        coords = new_coords;

        // Evaluate the new energy and force matrix
        new_energy = calcTotalEnergy();
        force_matrix = calcAnalyticalForces();

        // Change step size if norm of force matrix is too large
        if (old_energy < new_energy) {
            stepsize *= 0.5;
            coords = old_coords;
        }
        else {
            old_energy = new_energy;
            stepsize *= 1.1;
        }

        // Increment the iteration counter
        iter++;
    }

    // Print final geometry, energy, and norm of the force matrix
    cout << "Final Geometry after Steepest Descent:" << endl;
    print_geometry();
    cout << "Final LJ Energy: " << new_energy << endl;
    cout << "Final Norm of Forces: " << norm(force_matrix, "fro") << endl;
    cout << "Number of Iterations: " << iter << endl;

    // Return the final coordinates
    return coords;
}


    