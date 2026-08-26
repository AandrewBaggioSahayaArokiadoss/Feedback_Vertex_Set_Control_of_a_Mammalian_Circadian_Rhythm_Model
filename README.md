# Feedback Vertex Set Control of a 21-State Mammalian Circadian Rhythm Model

## Overview

This repository contains a Julia implementation of **feedback vertex set (FVS) control** for a nonlinear 21-state mammalian circadian rhythm model.

The model is represented by a directed state-interaction graph in which an edge from state `u` to state `v` means that `u` directly influences the differential equation governing `v`.

The script verifies six candidate feedback vertex sets, clamps the states in a selected FVS to their desired attractor values, simulates the remaining nonlinear dynamics, and checks whether the full system approaches the target long-term behavior.

## Main idea

For a chosen feedback vertex set, the corresponding state variables are initialized at their desired attractor values and kept fixed there during the simulation. The remaining states evolve according to the original nonlinear model.

The goal is to demonstrate that controlling only the FVS states is sufficient to drive the rest of the system toward the desired attractor.

The supplied target attractor is the 21-dimensional state vector

```text
[0.1229, 0.4361, 0.1040, 0.5496, 2.1465, 1.4523, 0.1700,
 5.0701, 0.7137, 0.3011, 2.6030, 1.0467, 1.4131, 0.2738,
 0.0971, 0.2326, 0.6503, 0.3251, 1.2300, 0.7011, 1.6326]
```

## State variables

The model contains the following 21 states:

```text
PER1, PER2, per1, per2, CRY1, CRY2, cry1, cry2,
PER1CRY1, PER1CRY2, PER2CRY1, PER2CRY2,
REVERB, reverb, RORC, rorc, CLK, clk, BMAL1, bmal1,
CLKBMAL1
```

The nonlinear model includes protein production and degradation, mRNA dynamics, protein-complex formation, and Hill-type activation and inhibition terms.

## Feedback vertex sets

The script contains six candidate FVSs:

1. `PER1, PER2, CRY1, CRY2, CLK, BMAL1, RORC`
2. `PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, RORC`
3. `PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, RORC`
4. `PER1, PER2, CRY1, CRY2, CLK, BMAL1, rorc`
5. `PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, rorc`
6. `PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, rorc`

Each candidate set is checked automatically by removing those vertices from the interaction digraph and verifying that the remaining graph is acyclic using a topological-sort procedure.

## Numerical simulation

The model is integrated with a fixed-step fourth-order Runge-Kutta method.

Default settings are:

```text
Final time: 25.0
Time step:  0.001
Initial state: all zeros, except controlled FVS states
Tolerance for convergence test: 1%
```

For each FVS, the code reports:

- whether the candidate is a valid feedback vertex set;
- convergence time under the 1% per-state relative-error criterion;
- final maximum relative error.

A CSV summary is also written to disk.

## Euclidean distance to the target trajectory

For the representative simulation, the code also computes the Euclidean distance between the full system state and the desired attractor:

```math
d(t) = \|x(t)-x^*\|_2
     = \sqrt{\sum_{i=1}^{21}(x_i(t)-x_i^*)^2}.
```

Here:

- `x(t)` is the 21-dimensional controlled system trajectory;
- `x*` is the desired attractor.

A decreasing value of `d(t)` indicates that the complete system state is moving toward the target attractor. If `d(t)` tends to zero, the controlled trajectory approaches the desired long-term state in the full 21-dimensional state space.

## Plots

The script uses **FVS 2** as the representative case for visualization and generates two figures.

### 1. All 21 state trajectories

```text
circadian_fvs2_all_state_trajectories.png
```

This figure contains 21 subplots, one for each state variable. In every subplot:

- the solid curve is the simulated state trajectory;
- the dashed horizontal line is the desired attractor value.

The controlled FVS states remain at their target values, while the remaining states evolve toward the attractor.

### 2. Euclidean distance to the attractor

```text
circadian_fvs2_euclidean_distance.png
```

This figure plots

```math
\|x(t)-x^*\|_2
```

against time and provides a single scalar measure of how close the complete 21-state system is to the desired attractor.

## Output files

Running the script generates:

```text
circadian_fvs_summary.csv
circadian_fvs2_all_state_trajectories.png
circadian_fvs2_euclidean_distance.png
```

## Running in Google Colab

The code is designed to run directly in a Julia-enabled Google Colab environment.

Upload the Julia file and run:

```julia
include("circadian_fvs_21states_euclidean_error.jl")
```

If `Plots.jl` is not already installed in the Julia environment, the script installs it automatically.

## Dependencies

- Julia
- `Plots.jl`
- `Printf` from the Julia standard library

## Structure of the code

The main module provides:

- `circadian_rhs` — nonlinear 21-state model;
- `circadian_controlled_rhs` — model with selected FVS derivatives clamped to zero;
- `rk4_simulate` — fourth-order Runge-Kutta simulation;
- `verify_fvs` — graph-based verification of a feedback vertex set;
- `convergence_time` — 1% relative-error convergence test;
- `max_relative_error` — maximum componentwise relative error;
- `euclidean_distance` — Euclidean distance of the full state vector from the desired attractor.

## Interpretation

The experiment illustrates the feedback-vertex-set control principle on a biologically motivated nonlinear dynamical system: instead of directly controlling all 21 states, only a structurally selected subset of states is constrained to follow the desired long-term behavior. The response of the remaining states is then observed to determine whether the full system converges toward the same attractor.
