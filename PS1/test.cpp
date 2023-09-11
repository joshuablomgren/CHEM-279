#include "gold_lj.h"

int main(){
    // Create a GoldCluster object and print the geometry
    GoldCluster cluster("example.txt");

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
    cout << "LJ energy of the cluster: " << total_energy << endl;

    // Test Analytical Force Calculation
    mat forces_analytical = cluster.calcAnalyticalForces();
    cout << "Analytical Forces: " << endl;
    cout << forces_analytical << endl;

    // Test Forward Difference Force Calculation
    mat forces_forward = cluster.forwardDifference(0.01);
    cout << "Forward Difference Forces (h = 0.01): " << endl;
    cout << forces_forward << endl;

    // Test Central Difference Force Calculation
    mat forces_central = cluster.centralDifference(0.01);
    cout << "Central Difference Forces (h = 0.01): " << endl;
    cout << forces_central << endl;

    // Separate the output
    cout << "----------------------------------------" << endl;

    // Test steepest descent
    mat coords_min = cluster.steepestDescent(0.01, 1e-10);

    return 0;
}