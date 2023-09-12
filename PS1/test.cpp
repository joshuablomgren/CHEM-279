#include "gold_lj.h"

int main(){
    cout << "1. Lennard Jones Energy" << endl;
    // Create a GoldCluster object and print the geometry
    GoldCluster cluster("4atoms.txt");

    cout << "Input Geometry: " << endl;
    cluster.print_geometry();

    // If non-gold atoms are present or file format is incorrect, return error
    for (int i = 0; i < cluster.num_atoms; i++) {
        if (cluster.z_vals[i] != 79) {
            cout << "Error: Non-gold atoms present in file!" << endl;
            return -1;
        }
        if (cluster.z_vals[i] <= 0) {
            cout << "Error: Invalid atomic number!" << endl;
            return -1;
        }
        if (cluster.coords.n_cols != 3) {
            cout << "Error: Coordinates are not 3D!" << endl;
            return -1;
        }
    }

    // Test Total Energy Calculation
    double total_energy = cluster.calcTotalEnergy();
    cout << "LJ energy of the cluster: " << total_energy << endl << endl;

    // Test again on 5 atom cluster
    GoldCluster cluster2("5atoms.txt");
    cout << "Input Geometry: " << endl;
    cluster2.print_geometry();
    double total_energy2 = cluster2.calcTotalEnergy();
    cout << "LJ energy of the cluster: " << total_energy2 << endl;

    cout << "----------------------------------------" << endl;
    cout << "2. Finite Difference vs Analytical Forces" << endl;

    cout << "Input Geometry: " << endl;
    cluster.print_geometry();

    // Test Analytical Force Calculation
    mat forces_analytical = cluster.calcAnalyticalForces();
    cout << endl << "Analytical Forces: " << endl;
    cout << forces_analytical << endl;

    // Test Forward Difference Force Calculation
    mat forces_forward = cluster.forwardDifference(0.01);
    cout << "Forward Difference Forces (h = 0.01): " << endl;
    cout << forces_forward << endl;

    // Test Central Difference Force Calculation
    mat forces_central = cluster.centralDifference(0.01);
    cout << "Central Difference Forces (h = 0.01): " << endl;
    cout << forces_central << endl;

    // Calculate Error for Forward Difference at different h values
    vec h_vals = {0.1, 0.01, 0.001, 0.0001};
    vec error_forward = zeros<vec>(4);
    for (int i = 0; i < 4; i++) {
        mat forces_forward = cluster.forwardDifference(h_vals(i));
        error_forward(i) = norm(forces_forward - forces_analytical);
    }

    // Calculate Error for Central Difference at different h values
    vec error_central = zeros<vec>(4);
    for (int i = 0; i < 4; i++) {
        mat forces_central = cluster.centralDifference(h_vals(i));
        error_central(i) = norm(forces_central - forces_analytical);
    }

    // Print the errors
    ofstream forward_file("forward_error.txt");
    cout << "Error for Forward Difference: " << endl;
    for (int i = 0; i < 4; i++) {
        cout << h_vals(i) << ": " << error_forward(i) << endl;
        forward_file << h_vals(i) << " " << error_forward(i) << endl;
    }
    forward_file.close();

    cout << endl;

    ofstream central_file("central_error.txt");
    cout << "Error for Central Difference: " << endl;
    for (int i = 0; i < 4; i++) {
        cout << h_vals(i) << ": " << error_central(i) << endl;
        central_file << h_vals(i) << " " << error_central(i) << endl;
    }
    central_file.close();

    // Separate the output
    cout << "----------------------------------------" << endl;
    cout << "3. Optimizing the Geometry of Small Gold Cluster" << endl;

    // Test steepest descent on 2 gold atoms
    GoldCluster cluster3("2atoms.txt");
    cluster3.steepestDescent(0.01, 1e-10);

    // Test steepest descent on 4 gold atoms
    cluster.steepestDescent(0.01, 1e-10);

    return 0;
}   