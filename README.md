# Feedback Vertex Set Control of a Mammalian Circadian Rhythm Model

**Short description:** Julia/Colab demonstration of feedback-vertex-set control in a 21-state mammalian circadian rhythm network, showing how clamping only FVS states to a desired natural attractor drives the remaining states toward that attractor.

## Overview

This repository reproduces, in Julia, a feedback-vertex-set (FVS) control experiment originally implemented in MATLAB for a nonlinear mammalian circadian rhythm model with 21 states.

The state-interaction digraph records which state variables influence the differential equations of which other state variables. For each candidate FVS, the code verifies that removing the selected states makes the remaining directed graph acyclic. The controlled FVS states are then set to their values at the desired natural attractor and held there, while all non-FVS states continue to evolve according to the original nonlinear model.

The numerical experiment tests whether controlling only the FVS coordinates is sufficient to drive the remaining states toward the desired long-term behavior.

## State ordering

```text
PER1, PER2, per1, per2, CRY1, CRY2, cry1, cry2,
PER1CRY1, PER1CRY2, PER2CRY1, PER2CRY2,
REVERB, reverb, RORC, rorc, CLK, clk, BMAL1, bmal1, CLKBMAL1
```

The desired attractor is

```text
[0.1229, 0.4361, 0.1040, 0.5496, 2.1465, 1.4523, 0.1700,
 5.0701, 0.7137, 0.3011, 2.6030, 1.0467, 1.4131, 0.2738,
 0.0971, 0.2326, 0.6503, 0.3251, 1.2300, 0.7011, 1.6326]
```

## Feedback vertex sets

The six FVS candidates used in the experiment are

```text
FVS 1 = PER1, PER2, CRY1, CRY2, CLK, BMAL1, RORC
FVS 2 = PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, RORC
FVS 3 = PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, RORC
FVS 4 = PER1, PER2, CRY1, CRY2, CLK, BMAL1, rorc
FVS 5 = PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, rorc
FVS 6 = PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, rorc
```

For every candidate, the code deletes the corresponding vertices and performs a topological-sort test. If all remaining vertices are included in the resulting topological order, the residual digraph is acyclic and the candidate is confirmed to be a feedback vertex set.

## Control experiment

The uncontrolled initial state is

```julia
x0 = zeros(21)
```

For a chosen FVS `S`, its coordinates are first placed at the desired attractor values:

```text
x_i(0) = x_i*,    i in S.
```

During the simulation those coordinates are clamped by imposing

```text
dx_i/dt = 0,      i in S,
```

while every non-FVS coordinate continues to follow its original differential equation.

The default numerical settings are

```text
final time = 25
time step  = 0.001
tolerance  = 1%
```

A fixed-step fourth-order Runge-Kutta scheme is used.

## What the code produces

Running the file:

1. verifies all six proposed FVSs;
2. simulates FVS control for all six choices;
3. reports the convergence time and final maximum relative error numerically;
4. writes `circadian_fvs_summary.csv`;
5. plots the normalized non-FVS trajectories for FVS 2;
6. plots all 21 individual state trajectories for FVS 2.

The generated figures are

```text
circadian_fvs2_normalized_nonfvs_states.png
circadian_fvs2_all_state_trajectories.png
```

The first figure divides every uncontrolled state by its desired attractor value, so convergence corresponds to each curve approaching `1`. The second figure shows each state in its original units together with its target value as a dashed horizontal line.

The previous maximum-relative-error figure is intentionally not generated. Relative errors are still calculated internally because they are useful for the convergence test and numerical summary.

## Repository structure

```text
mammalian-circadian-fvs-julia/
├── README.md
├── Project.toml
└── circadian_fvs_colab.jl
```

The Julia file is intentionally self-contained so that it can be uploaded directly to a Julia-enabled Google Colab session.

## Run in Google Colab

If Julia is already available in your Colab environment, upload `circadian_fvs_colab.jl` and execute

```julia
include("circadian_fvs_colab.jl")
```

The script assumes Julia itself is present. If `Plots.jl` is missing from the active Julia environment, it installs that package automatically.

## Run locally

From the repository directory:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. circadian_fvs_colab.jl
```

## MATLAB-to-Julia correspondence

| MATLAB file/function | Julia counterpart |
|---|---|
| `mammcircfunc.m` | `rhs` |
| `mammcircoverride.m` | `controlled_rhs` and exact FVS clamp |
| `mammcirc_graph_init.m` | `INFLUENCERS` |
| `mammcirc_fvs_init.m` | `FVS_SETS` |
| `name_position.m` | `state_indices` |
| `tolconverge.m` | `convergence_time` |
| `mammcirccheck.m` | experiment section of `circadian_fvs_colab.jl` |

## Scope

This repository is a computational reproduction of the supplied FVS-control experiment. The theoretical assumptions of the feedback-vertex-set control theorem are not re-derived here. The purpose of the code is to demonstrate numerically, for this 21-state circadian model, how prescribing only a structurally selected subset of state variables can drive the remaining network toward a desired natural attractor.
