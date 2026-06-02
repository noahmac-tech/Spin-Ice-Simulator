#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <future>
#include <thread>
#include <algorithm>

using namespace std;

// Struct to hold individual spin data
struct Spin {
    int pointing;
    double x, y, z;
    double X, Y, Z;
    int neighbours[6];
};

// Struct to hold the results returned by each thread
struct SimulationResult {
    double T0;
    double average_E;
    double average_M;
    double specific_heat;
    double susceptibility;
};

// SpinIceSimulation Class
class SpinIceSimulation {
private:
    int L;
    int num_spins;
    double a;
    double J_eff, Bx, By, Bz;
    vector<Spin> lattice;

    static constexpr double basis[16][3] = {
        {0.50, 0.50, 0.50}, {0.50, 0.00, 0.00}, {0.00, 0.50, 0.00}, {0.00, 0.00, 0.50},
        {0.50, 0.75, 0.75}, {0.00, 0.75, 0.25}, {0.50, 0.25, 0.25}, {0.00, 0.25, 0.75},
        {0.25, 0.75, 0.00}, {0.25, 0.25, 0.50}, {0.75, 0.75, 0.50}, {0.75, 0.25, 0.00},
        {0.25, 0.50, 0.25}, {0.75, 0.50, 0.75}, {0.25, 0.00, 0.75}, {0.75, 0.00, 0.25}
    };

    static constexpr double inv_sqrt3 = 0.57735026919;
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

                if (abs((dx * dx + dy * dy + dz * dz) - nn_dist_sq) < tolerance) {
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
        for (int i = 0; i < num_spins; ++i) { total_M_z += lattice[i].pointing * lattice[i].Z; }
        return total_M_z / num_spins;
    }

    double calculate_energy() const {
        double total_E = 0;
        for (int i = 0; i < num_spins; ++i) {
            double local_exchange = 0;
            for (int n = 0; n < 6; ++n) { local_exchange += lattice[lattice[i].neighbours[n]].pointing; }
            total_E += ((J_eff / 3.0) * lattice[i].pointing * local_exchange) / 2.0;
            double dot_product = (Bx * lattice[i].X) + (By * lattice[i].Y) + (Bz * lattice[i].Z);
            total_E -= lattice[i].pointing * dot_product;
        }
        return total_E / num_spins;
    }

    void metropolis_sweep(double beta, mt19937& gen, uniform_int_distribution<>& random_spin, uniform_real_distribution<>& probability) {
        for (int step = 0; step < num_spins; ++step) {
            int i = random_spin(gen); 
            int sum_neighbours = 0;
            for (int n = 0; n < 6; ++n) { sum_neighbours += lattice[lattice[i].neighbours[n]].pointing; }
            double delta_E_exchange = -2.0 * (J_eff / 3.0) * lattice[i].pointing * sum_neighbours;
            double dot_product = (Bx * lattice[i].X) + (By * lattice[i].Y) + (Bz * lattice[i].Z);
            double delta_E_zeeman = 2.0 * lattice[i].pointing * dot_product; 
            double delta_E = delta_E_exchange + delta_E_zeeman;

            if (delta_E < 0 || probability(gen) < exp(-delta_E * beta)) {
                lattice[i].pointing *= -1;
            }
        }
    }
};

SimulationResult run_simulation_for_temperature(double T0, int L, double a, double J_eff, double Bx, double By, double Bz, int eq_sweeps, int meas_sweeps, int sweeps_between) {
    
    // CRITICAL: Each thread needs its own random number generator!
    // We seed it using the current time + the temperature so they are all unique.
    random_device rd;
    mt19937 gen(rd() + static_cast<unsigned int>(T0 * 1000)); 
    uniform_real_distribution<> probability(0.0, 1.0);

    // Initialize the simulation object specifically for this thread
    SpinIceSimulation sim(L, a, J_eff, Bx, By, Bz);
    sim.build_lattice(gen);
    uniform_int_distribution<> random_spin(0, sim.get_num_spins() - 1);

    double beta = 1.0 / T0;

    // Thermalization
    for (int sweep = 0; sweep < eq_sweeps; ++sweep) {
        sim.metropolis_sweep(beta, gen, random_spin, probability);
    }

    // Measurement
    double magnetisation_sum = 0, magnetisation_sum_squared = 0;
    double energy_sum = 0, energy_sum_squared = 0;

    for (int measurement = 0; measurement < meas_sweeps; ++measurement) {
        for (int sweep = 0; sweep < sweeps_between; ++sweep){
            sim.metropolis_sweep(beta, gen, random_spin, probability);
        }
        double M = sim.calculate_magnetisation();
        double E = sim.calculate_energy();
        magnetisation_sum += M;
        magnetisation_sum_squared += M * M;
        energy_sum += E;
        energy_sum_squared += E * E;
    }

    // Compute averages
    double average_E = energy_sum / meas_sweeps;
    double average_E_squared = energy_sum_squared / meas_sweeps;
    double var_E = average_E_squared - (average_E * average_E);
    double specific_heat = (sim.get_num_spins() / (T0 * T0)) * var_E;
    
    double average_M = magnetisation_sum / meas_sweeps;
    double average_M_squared = magnetisation_sum_squared / meas_sweeps;
    double var_M = average_M_squared - (average_M * average_M);
    double susceptibility =  (sim.get_num_spins() / T0) * var_M;

    cout << "Completed T = " << T0 << " on thread ID: " << this_thread::get_id() << "\n";

    return {T0, average_E, average_M, specific_heat, susceptibility};
}

int main() {
    const int L = 4;
    const double a = 19.22495270;
    const double J_eff = 1.11;
    const double Bx = 0.0, By = 0.0, Bz = 2.0;
    const int eq_sweeps = 50000;
    const int meas_sweeps = 5000;
    const int sweeps_between = 20;

    // A vector to hold our "Futures" (promises that a thread will eventually return a SimulationResult)
    vector<future<SimulationResult>> futures;

    cout << "Launching multithreaded Spin Ice simulations across all CPU cores...\n";

    // Dispatch a thread for every single temperature
    for (double T0 = 6.0; T0 >= 0.1; T0 -= 0.1) {
        futures.push_back(
            async(launch::async, run_simulation_for_temperature, T0, L, a, J_eff, Bx, By, Bz, eq_sweeps, meas_sweeps, sweeps_between)
        );
    }

    // Wait for all threads to finish and collect the results
    vector<SimulationResult> results;
    for (auto& f : futures) {
        results.push_back(f.get()); // .get() blocks until the specific thread is done
    }

    // Sort the results (Because threads finish in random order, the data needs sorting!)
    sort(results.begin(), results.end(), [](const SimulationResult& a, const SimulationResult& b) {
        return a.T0 < b.T0; // Sort from lowest to highest temperature
    });

    // Write sorted results to file
    fstream all_data_file;
    all_data_file.open("data/data_spin_ice_parallel.txt", ios::out);
    all_data_file << "T0\tAverage_Energy\tAverage_Magnetisation\tSpecific_Heat\tSusceptibility\n";

    for (const auto& res : results) {
        all_data_file << res.T0 << "\t" << res.average_E << "\t" << res.average_M << "\t" 
                      << res.specific_heat << "\t" << res.susceptibility << "\n";
    }
    
    all_data_file.close();
    cout << "All threads finished successfully! Data saved to data/data_spin_ice_parallel.txt\n";
    return 0;
}