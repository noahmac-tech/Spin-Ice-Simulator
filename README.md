# Spin Ice & 2D Ising Model Monte Carlo Simulations

[![C++](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-yellow.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains high-performance Monte Carlo simulations of statistical mechanics models, specifically focusing on the **3D Spin Ice model**. 

The core physics engine is written in Object-Oriented C++ to leverage high-speed computation for millions of Metropolis-Hastings sweeps. The resulting thermodynamic data is processed, analyzed, and visualized using a custom Python command-line interface (CLI).

## Key Features
* **Metropolis-Hastings Algorithm:** Robust implementation of Markov Chain Monte Carlo (MCMC) methods to simulate spin dynamics, thermalization, and phase transitions.
* **3D Spin Ice Thermodynamics:** Generates a 16-basis tetrahedral lattice to calculate Energy, Magnetisation, Specific Heat ($C_v$), and Magnetic Susceptibility ($\chi$) across varying temperatures.
* **Residual Entropy Calculation:** Simulates the highly degenerate "two-in, two-out" ground state of Spin Ice at near-zero Kelvin (in zero magnetic field) to observe Pauling's residual entropy plateau.
* **Data Visualization CLI:** A Python `argparse` tool utilizing Pandas and Matplotlib to cleanly generate publication-ready plots.

---

## Project Structure

```text
Spin-Ice-Simulation/
├── src/
│   ├── spin_ice.cpp             # Main thermodynamics simulation
│   └── low_temp_entropy.cpp     # Zero B-field entropy simulation
├── scripts/
│   └── plot_spin_ice.py         # Data visualization CLI
├── data/                        # Output datasets (.txt)
├── figures/                     # Generated plots (.png)
├── requirements.txt             # Python dependencies
└── README.md