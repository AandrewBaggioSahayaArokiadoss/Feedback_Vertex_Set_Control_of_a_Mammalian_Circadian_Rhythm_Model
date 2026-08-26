# README: Mammalian Circadian Rhythm FVS Control

## Overview
This project implements and verifies a Feedback Vertex Set (FVS) control strategy for a 21-state ODE model of the mammalian circadian rhythm. The goal is to demonstrate that by controlling a specific subset of variables (the FVS), the entire system can be driven toward a desired steady-state attractor.

## Features
- **21-State Circadian Model**: A Julia translation of a complex mammalian circadian rhythm model involving PER, CRY, CLK, BMAL1, and other regulatory proteins/complexes.
- **Influence Digraph Analysis**: Construction of the interaction network to identify dependencies between states.
- **FVS Verification**: Implementation of Kahn's algorithm to verify that proposed sets are indeed Feedback Vertex Sets (i.e., their removal leaves the remaining graph acyclic).
- **RK4 Simulation**: A 4th-order Runge-Kutta integrator with a clamping mechanism to simulate controlled trajectories.
- **Convergence Analysis**: Tools to calculate relative error and identify the exact time at which the system converges to the attractor within a 1% tolerance.

## Control Strategy
The control is implemented by "clamping" the FVS states. In this simulation, the states identified in the FVS are initialized at their target attractor values and held constant ($\\dot{x}_i = 0$), forcing the remaining variables to follow the acyclic flow toward the equilibrium.

## Results
The simulation tests six different FVS candidates. Key metrics include:
- **Convergence Time**: The time required for all 21 states to reach the attractor.
- **Max Relative Error**: The final precision of the simulation compared to the target steady state.

## Requirements
- Julia 1.x
- `Plots.jl` for visualization
- `Printf` for formatted output
