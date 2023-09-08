#include "gold_lj.h"

int main(){
    // Create a GoldCluster object and print the geometry
    GoldCluster cluster("example.txt");
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

    // Test Distance Calculation
    double distance = cluster.calcDistance(0, 1);
    cout << "Distance between atoms 0 and 1: " << distance << endl;

    // Test LJ Calculation
    double lj = cluster.calcLJ(0, 1);
    cout << "LJ energy between atoms 0 and 1: " << lj << endl;

    // Test Total Energy Calculation
    double total_energy = cluster.calcTotalEnergy();
    cout << "Total energy of the cluster: " << total_energy << endl;

    return 0;
}