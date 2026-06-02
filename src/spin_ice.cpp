#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

using namespace std;

// Struct to hold individual spin data
struct Spin {
    int pointing;
    double x, y, z;
    double X, Y, Z;
    int neighbours[6];
};

// SpinIceSimulation Class encapsulates the data and the physics
class SpinIceSimulation {
private:
    int L;
    int num_spins;
    double a;
    double J_eff, Bx, By, Bz;
    vector<Spin> lattice;

    // Static constexpr arrays for the lattice basis and axes
    static constexpr double basis[16][3] = {
        {0.50, 0.50, 0.50}, {0.50, 0.00, 0.00}, {0.00, 0.50, 0.00}, {0.00, 0.00, 0.50},
        {0.50, 0.75, 0.75}, {0.00, 0.75, 0.25}, {0.50, 0.25, 0.25}, {0.00, 0.25, 0.75},
        {0.25, 0.75, 0.00}, {0.25, 0.25, 0.50}, {0.75, 0.75, 0.50}, {0.75, 0.25, 0.00},
        {0.25, 0.50, 0.25}, {0.75, 0.50, 0.75}, {0.25, 0.00, 0.75}, {0.75, 0.00, 0.25}
    };

    static constexpr double inv_sqrt3 = 0.57735026919; // 1.0 / sqrt(3.0)
    static constexpr double axes[16][3] = {
        {inv_sqrt3, inv_sqrt3, inv_sqrt3},   {inv_sqrt3, -inv_sqrt3, -inv_sqrt3}, 
        {-inv_sqrt3, inv_sqrt3, -inv_sqrt3}, {-inv_sqrt3, -inv_sqrt3, inv_sqrt3},
        {inv_sqrt3, -inv_sqrt3, -inv_sqrt3}, {-inv_sqrt3, -inv_sqrt3, inv_sqrt3}, 
        {inv_sqrt3, inv_sqrt3, inv_sqrt3},   {-inv_sqrt3, inv_sqrt3, -inv_sqrt3},
        {-inv_sqrt3, inv_sqrt3, -inv_sqrt3}, {inv_sqrt3, inv_sqrt3, inv_sqrt3},   
        {-inv_sqrt3, -inv_sqrt3, inv_sqrt3}, {inv_sqrt3, -inv_sqrt3, -inv_sqrt3},
        {-inv_sqrt3, -inv_sqrt3, inv_sqrt3}, {-inv_sqrt3, inv_sqrt3, -inv_sqrt3}, 
        {inv_sqrt3, -inv_sqrt3, -inv_sqrt3}, {inv_sqrt3, inv_sqrt3, inv_sqrt3}
    };

public:
    // Constructor initializes the parameters and allocates the lattice
    SpinIceSimulation(int size, double lattice_constant, double J, double bx, double by, double bz) 
        : L(size), a(lattice_constant), J_eff(J), Bx(bx), By(by), Bz(bz) {
        num_spins = L * L * L * 16;
        lattice.resize(num_spins);
    }

    int get_num_spins() const { return num_spins; }

    void build_lattice(mt19937& gen) {
        uniform_int_distribution<> random_num(0, 1);
        int counter = 0;
        double width = L * a;
        
        // Assign positions and spin axes
        for (int ix = 0; ix < L; ++ix) {
            for (int iy = 0; iy < L; ++iy) {
                for (int iz = 0; iz < L; ++iz) {
                    for (int s = 0; s < 16; ++s) {
                        lattice[counter].pointing = (random_num(gen) == 0) ? -1 : 1;
                        lattice[counter].x = (ix + basis[s][0]) * a;
                        lattice[counter].y = (iy + basis[s][1]) * a;
                        lattice[counter].z = (iz + basis[s][2]) * a;
                        lattice[counter].X = axes[s][0];
                        lattice[counter].Y = axes[s][1];
                        lattice[counter].Z = axes[s][2];
                        counter++;
                    }
                }
            }
        }

        // Map Nearest Neighbours using periodic boundary conditions
        double nn_dist_sq = 0.125 * a * a;
        double tolerance = 1e-4 * a * a; 
        
        for (int i = 0; i < num_spins; ++i) {
            int neighbour_count = 0;
            for (int j = 0; j < num_spins; ++j) {
                if (i == j) continue;

                double dx = lattice[i].x - lattice[j].x;
                double dy = lattice[i].y - lattice[j].y;
                double dz = lattice[i].z - lattice[j].z;

                dx -= width * round(dx / width);
                dy -= width * round(dy / width);
                dz -= width * round(dz / width);

                double distance_squared = dx * dx + dy * dy + dz * dz;

                if (abs(distance_squared - nn_dist_sq) < tolerance) {
                    if (neighbour_count < 6) {
                        lattice[i].neighbours[neighbour_count] = j;
                        neighbour_count++;
                    }
                }
            }
        }
    }

    double calculate_magnetisation() const {
        double total_M_z = 0;
        for (int i = 0; i < num_spins; ++i) {
            total_M_z += lattice[i].pointing * lattice[i].Z; 
        }
        return total_M_z / num_spins;
    }

    double calculate_energy() const {
        double total_E = 0;
        for (int i = 0; i < num_spins; ++i) {
            double local_exchange = 0;
            for (int n = 0; n < 6; ++n) {
                int neighbour = lattice[i].neighbours[n];
                local_exchange += lattice[neighbour].pointing;
            }
            total_E += ((J_eff / 3.0) * lattice[i].pointing * local_exchange) / 2.0;

            double dot_product = (Bx * lattice[i].X) + (By * lattice[i].Y) + (Bz * lattice[i].Z);
            total_E -= lattice[i].pointing * dot_product;
        }
        return total_E / num_spins;
    }

    // Sweep Logic for Metropolis Algorithm
    void metropolis_sweep(double beta, mt19937& gen, uniform_int_distribution<>& random_spin, uniform_real_distribution<>& probability) {
        for (int step = 0; step < num_spins; ++step) {
            int i = random_spin(gen); 
            int sum_neighbours = 0;
            for (int n = 0; n < 6; ++n) {
                sum_neighbours += lattice[lattice[i].neighbours[n]].pointing;
            }
            double delta_E_exchange = -2.0 * (J_eff / 3.0) * lattice[i].pointing * sum_neighbours;
            double dot_product = (Bx * lattice[i].X) + (By * lattice[i].Y) + (Bz * lattice[i].Z);
            double delta_E_zeeman = 2.0 * lattice[i].pointing * dot_product; 
            double delta_E = delta_E_exchange + delta_E_zeeman;

            if (delta_E < 0 || probability(gen) < exp(-delta_E * beta)) {
                lattice[i].pointing *= -1; // Flip accepted
            }
        }
    }
};

int main() {
    // Configuration Parameters
    const int L = 4;
    const double a = 19.22495270;
    const double J_eff = 1.11;
    const double Bx = 0.0, By = 0.0, Bz = 2.0;
    const int equilibration_sweeps = 50000;
    const int measurement_sweeps = 5000;
    const int sweeps_between_measurements = 20;

    // Setup RNG and initialize simulation object
    random_device rd;
    mt19937 gen(rd()); 
    uniform_real_distribution<> probability(0.0, 1.0);

    SpinIceSimulation sim(L, a, J_eff, Bx, By, Bz);
    sim.build_lattice(gen);

    uniform_int_distribution<> random_spin(0, sim.get_num_spins() - 1);

    // Setup output file
    fstream all_data_file;
    all_data_file.open("data/data_spin_ice.txt", ios::out);
    all_data_file << "T0\tAverage_Energy\tAverage_Magnetisation\tSpecific_Heat\tSusceptibility\n";

    // Run Simulation over Temperatures
    for (double T0 = 6.0; T0 >= 0.1; T0 -= 0.1) {
        double beta = 1.0 / T0;
        cout << "Running T = " << T0 << "..." << flush;

        // Thermalization / Equilibration Phase
        for (int sweep = 0; sweep < equilibration_sweeps; ++sweep) {
            sim.metropolis_sweep(beta, gen, random_spin, probability);
        }

        // Measurement Phase
        double magnetisation_sum = 0, magnetisation_sum_squared = 0;
        double energy_sum = 0, energy_sum_squared = 0;

        for (int measurement = 0; measurement < measurement_sweeps; ++measurement) {
            // Decorrelate before taking a measurement
            for (int sweep = 0; sweep < sweeps_between_measurements; ++sweep){
                sim.metropolis_sweep(beta, gen, random_spin, probability);
            }
            
            double M = sim.calculate_magnetisation();
            double E = sim.calculate_energy();

            magnetisation_sum += M;
            magnetisation_sum_squared += M * M;
            energy_sum += E;
            energy_sum_squared += E * E;
        }

        // Compute statistical averages
        double average_E = energy_sum / measurement_sweeps;
        double average_E_squared = energy_sum_squared / measurement_sweeps;
        double var_E = average_E_squared - (average_E * average_E);
        double specific_heat = (sim.get_num_spins() / (T0 * T0)) * var_E;
        
        double average_M = magnetisation_sum / measurement_sweeps;
        double average_M_squared = magnetisation_sum_squared / measurement_sweeps;
        double var_M = average_M_squared - (average_M * average_M);
        double susceptibility =  (sim.get_num_spins() / T0) * var_M;

        all_data_file << T0 << "\t" << average_E << "\t" << average_M << "\t" << specific_heat << "\t" << susceptibility << "\n";
        cout << " Done.\n";
    }
    
    all_data_file.close();
    cout << "Script is finished!\n";
    return 0;
}