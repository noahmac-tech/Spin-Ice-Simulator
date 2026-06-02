import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# Set global plotting parameters
plt.rcParams.update({'font.size': 10})

# Ensure the output directory for figures exists
os.makedirs("figures", exist_ok=True)

def plot_spin_ice_thermodynamics(data_file="data/data_spin_ice.txt"):
    """
    Reads Spin Ice simulation data and plots Energy, Specific Heat, 
    Magnetisation, and Susceptibility as a function of Temperature.
    """
    try:
        data = pd.read_csv(data_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Could not find '{data_file}'. Have you run the simulation yet?")
        return

    data = data.sort_values(by='T0')
    
    # --- FIGURE 1: Energy and Specific Heat ---
    fig1, axes1 = plt.subplots(2, 1, figsize=(4.0, 6.0), sharex=True)
    
    axes1[0].plot(data['T0'], data['Average_Energy'], marker='o', color='blue', linestyle='None', markersize=4)
    axes1[0].set_ylabel(r'Energy per spin $\langle E/N \rangle$')
    axes1[0].grid(True, linestyle='--', alpha=0.6)
    
    axes1[1].plot(data['T0'], data['Specific_Heat'], marker='^', color='red', linestyle='None', markersize=4)
    axes1[1].set_xlabel('Temperature (K)')
    axes1[1].set_ylabel('Specific Heat ($c / k_B$)')
    axes1[1].grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    fig1_path = "figures/Spin_Ice_Energy_Cv.png"
    plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {fig1_path}")
    
    # --- FIGURE 2: Magnetisation and Susceptibility ---
    fig2, axes2 = plt.subplots(2, 1, figsize=(4.0, 6.0), sharex=True)
    
    axes2[0].plot(data['T0'], data['Average_Magnetisation'], marker='s', color='green', linestyle='None', markersize=4)
    axes2[0].set_ylabel('Magnetisation ($M_z$)')
    axes2[0].grid(True, linestyle='--', alpha=0.6)
    
    axes2[1].plot(data['T0'], data['Susceptibility'], marker='d', color='purple', linestyle='None', markersize=4)
    axes2[1].set_xlabel('Temperature (K)')
    axes2[1].set_ylabel(r'Susceptibility ($\chi$)')
    axes2[1].grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    fig2_path = "figures/Spin_Ice_Mag_Sus.png"
    plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {fig2_path}")


def plot_low_temp_entropy(data_file="data/data_low_temp_entropy.txt"):
    """
    Calculates and plots the residual entropy of the Spin Ice model at low temperatures 
    by integrating the specific heat.
    """
    try:
        data = pd.read_csv(data_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Could not find '{data_file}'.")
        return

    data = data.sort_values(by='T0')
    
    T = data['T0']
    specific_heat = data['Specific_Heat']
    
    # Calculate entropy via numerical integration: S(T) = \int (C_v / T) dT
    entropy = cumulative_trapezoid(specific_heat / T, T, initial=0)
    
    fig, axes = plt.subplots(2, 1, figsize=(4.0, 6.0), sharex=True)
    
    # Top Plot: Specific Heat
    axes[0].plot(T, specific_heat, marker='^', color='red', linestyle='None', markersize=4)
    axes[0].set_ylabel('Specific Heat ($c / k_B$)')
    axes[0].grid(True, linestyle='--', alpha=0.6)
    
    # Bottom Plot: Integrated Entropy
    axes[1].plot(T, entropy, marker='o', color='purple', linestyle='None', markersize=4, label='Simulation')
    axes[1].set_xlabel('Temperature (K)')
    axes[1].set_ylabel(r'$\Delta S(T)$ per spin')
    axes[1].grid(True, linestyle='--', alpha=0.6)
        
    # Pauling's residual entropy estimation
    max_entropy = np.log(2)
    residual_entropy = 0.5 * np.log(1.5)
    plateau = max_entropy - residual_entropy
    
    axes[1].axhline(y=max_entropy, color='black', linestyle=':', label=r'Total possible entropy: $\ln(2)$')
    axes[1].axhline(y=plateau, color='green', linestyle='--', label=f'Expected Plateau: {plateau:.3f}')
    axes[1].legend(loc='lower right', fontsize=8)

    plt.tight_layout()
    fig_path = "figures/Spin_Ice_Entropy.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate plots for the Spin Ice Monte Carlo Simulation.")
    parser.add_argument('--plot', type=str, choices=['thermo', 'entropy', 'all'], default='all',
                        help="Choose which plot to generate: 'thermo', 'entropy', or 'all'.")
    parser.add_argument('--parallel', action='store_true', 
                        help="Read data from the multithreaded simulations (*_parallel.txt)")
    parser.add_argument('--show', action='store_true', 
                        help="Flag to display the plots interactively after saving.")

    args = parser.parse_args()

    print("Generating plots...")
    
    # Determine which files to read and what to name the output images based on the parallel flag
    if args.parallel:
        thermo_data = "data/data_spin_ice_parallel.txt"
        entropy_data = "data/data_low_temp_entropy_parallel.txt"
        prefix = "Spin_Ice_Parallel"
    else:
        thermo_data = "data/data_spin_ice.txt"
        entropy_data = "data/data_low_temp_entropy.txt"
        prefix = "Spin_Ice"

    # Call the plotting functions with dynamic variables
    if args.plot in ['thermo', 'all']:
        plot_spin_ice_thermodynamics(data_file=thermo_data, output_prefix=prefix)
        
    if args.plot in ['entropy', 'all']:
        plot_low_temp_entropy(data_file=entropy_data, output_prefix=prefix)

    if args.show:
        plt.show()