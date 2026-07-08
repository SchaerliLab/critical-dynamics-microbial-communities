# Critical Dynamics in Microbial Communities

This repository contains MATLAB scripts for simulating dynamical models of competing microbial populations coupled through slowly evolving environmental variables.

The model can exhibit different dynamical regimes depending on the interaction parameters:

- self-organized criticality (SOC),
- self-organized bistability (SOB),
- oscillatory dynamics.

---

## Repository Structure

Both scripts include deterministic and stochastic simulations and allow exploration of three dynamical regimes: self-organized criticality (SOC), self-organized bistability (SOB), and oscillations. Example parameter sets for each regime can be selected by setting

```matlab
behavior = 'oscillation'; % 'oscillation', 'SOC', 'SOB'
```

The scripts also generate the figures used to analyze the dynamics, including event-size distributions, bifurcation diagrams, and stochastic and deterministic trajectories.

### `ODE_simulation.m`

Implements the reduced **A-B-a model**, describing two competing populations (`A` and `B`) coupled through a slowly evolving environmental variable (`a`).

In addition to numerical simulations, the script computes the effective potential of the reduced system. Setting

```matlab
show_simplified = true;
```

visualizes the deterministic trajectories of both the original **A-B-a** model and its reduced **molar ratio-a** approximation.

### `ODE_simulation_full.m`

Implements the full **A-B-a-b model**, including two slowly evolving environmental variables (`a` and `b`).

The script computes the phase-space nullclines and overlays stochastic and deterministic trajectories. Setting

```matlab
show_trajectory = true;
```

visualizes the full deterministic trajectory of the system.

---

## Requirements

- MATLAB R2022b or newer
- MATLAB Symbolic Math Toolbox

---

## Citation

If you use this code, please cite:

G. Holló, J. H. Park, I. Lagzi, Y. Schaerli,  
*Critical dynamics shape stability in microbial communities*.
