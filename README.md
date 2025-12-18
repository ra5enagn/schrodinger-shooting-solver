# Schrödinger Equation Solver: Shooting Method

This repository contains a MATLAB implementation of the **Shooting Method**, a numerical technique used to find the bound-state eigenenergies and corresponding wavefunctions for a finite potential well.

## Overview

The Shooting Method converts the Time-Independent Schrödinger Equation (TISE), a boundary value problem, into an initial value problem. By integrating the wavefunction from both ends of a spatial mesh and comparing them at a matching point, we can iteratively determine the energies ($E$) that allow for a physically valid, continuous wavefunction.

$$-\frac{\hbar^2}{2m^*} \frac{d^2\psi}{dx^2} + V(x)\psi = E\psi$$

## Physics Configuration

The simulation is currently set up for a **GaAs (Gallium Arsenide)** quantum well with the following parameters:
* **Effective Mass ($m^*$):** $0.2 \cdot m_0$ (where $m_0$ is the free electron mass).
* **Well Width:** $10$ nm.
* **Well Depth:** $1$ eV.
* **Mesh Resolution:** $1500$ points total ($500$ points for each barrier and the well).



## How the Algorithm Works

1.  **Potential Initialization:** A spatial mesh is defined where the potential $V(x)$ is $1$ eV in the barriers and $0$ eV in the well.
2.  **Forward & Backward Integration:** For a given energy guess ($E_{guess}$), the code integrates the wavefunction $\psi$ starting from $\psi=0$ at the far left and far right boundaries.
3.  **Matching:** The two solutions are compared at a `matching_point`. The ratio between the two is used to normalize the backward "shot" to the forward "shot".
4.  **Error Calculation:** The discontinuity between the two wavefunctions at the point adjacent to the matching point is recorded as the "Matching Error".
5.  **Eigenvalue Search:** The code sweeps through energy levels. Local minima in the error curve indicate that a valid Eigen-energy has been found.
6.  **Reconstruction:** The script re-calculates the wavefunctions specifically for the detected eigenenergies for final visualization.

## Visualizations

The script produces two primary figures:
* **Figure 1:** The Potential Profile overlaid with the calculated wavefunctions, shifted and labeled by their specific Eigen-energies (eV).
* **Figure 2:** The Matching Error vs. Energy curve, showing the sharp dips that represent valid quantum states.



## Requirements
* MATLAB (R2018b or later recommended)
* No additional toolboxes required.
