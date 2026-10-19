# 3D Printer Hybrid System Simulation

This repository contains a Python-based simulation for a 3D printer modeled as a hybrid automaton. The project integrates continuous physical dynamics with discrete supervisory logic to accurately model mechanical motion, electrical circuits, and thermal properties during a print interruption and recovery cycle.

## System Architecture

The simulation is driven by a 10-dimensional continuous vector field and bounded by discrete event guards using `scipy.integrate.solve_ivp`.

* **Continuous Dynamics (10D State Space):** Models DC motor positions, velocities, and currents, alongside the nozzle's RLC circuit, thermal temperature, and an integrated odometer to track filament usage and spatial displacement.
* Finite State Machine:
  * **`q0`**: head is moving moving forward.
  * **`q1`**: save the current X and Y coordinates.
  * **`q2`**: head returning to the wipe tower at the origin.
  * **`q3`**: Scrubbing and wiping sequence at the origin.
  * **`q4`**: Thermal cooling phase targeting a safe threshold of 30°C.
  * **`q5`**: 5-second filament loading delay.
  * **`q6`**: Thermal heating phase targeting an operational threshold of 200°C.
  * **`q7`**: return to the saved X and Y print coordinates.

## Dependencies

The simulation requires standard scientific Python libraries:
* `numpy`
* `scipy`
* `matplotlib`

## Usage

Execute the main simulation script from the terminal:

```bash
python3 hybrid_system.py
```

The script will output timestamped discrete event triggers to the console and generate a comprehensive `matplotlib` figure containing the trajectory and state data.

## Simulation Results

![Hybrid System Simulation Plots] 
<img width="4200" height="2700" alt="simulation_hybrid" src="https://github.com/user-attachments/assets/3e4707b9-6817-42ee-a60e-3d793a8eb0b4" />


The generated output successfully demonstrates the complete hybrid cycle:
* **Cartesian Trajectory:** Shows the precise open-loop path tracking the tool head from the origin, stopping at the saved position (5.00, 0.00), reversing to the Wipe Tower, and returning to the saved position.
* **Nozzle Temperature:** Validates the cooling curve down to 30°C and the heating curve past 200°C.
* **FSM Transitions:** Confirms the strict chronological progression of the 8 discrete modes (q0 through q7) against the master simulation clock.
