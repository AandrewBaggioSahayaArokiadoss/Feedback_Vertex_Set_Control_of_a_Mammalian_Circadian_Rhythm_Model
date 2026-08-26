# Feedback Vertex Set Control of a Mammalian Circadian Rhythm Model

**Repository short description:** Julia/Colab reproduction of feedback-vertex-set control in a 21-state mammalian circadian rhythm network, showing how clamping only FVS states to a desired natural attractor drives the remaining states toward that attractor.

## What this repository demonstrates

The supplied MATLAB code implements a 21-state nonlinear mammalian circadian rhythm model together with an interaction digraph describing which states influence which other states. It also contains six candidate feedback vertex sets (FVSs).

The computational idea is:

1. Build the directed state-interaction graph.
2. Choose a feedback vertex set, i.e. a set of vertices whose removal makes the remaining directed graph acyclic.
3. Select a desired natural long-term behavior of the original system.
4. Make only the FVS states follow their values on that desired behavior.
5. Leave every non-FVS state uncontrolled.
6. Simulate the full nonlinear model and observe whether the remaining states are driven toward the desired behavior.

In the uploaded MATLAB experiment, the desired long-term behavior is a **fixed attractor/equilibrium vector**. Therefore the FVS states are initialized at their attractor values and then clamped there by setting their derivatives to zero. This repository reproduces that experiment in Julia.

> The control-theoretic theorem motivating feedback-vertex-set control has assumptions that are not contained in the uploaded MATLAB files. This repository is therefore a numerical reproduction of the supplied experiment, not a derivation or verification of the theorem itself.

## Model

The state ordering is

```text
PER1, PER2, per1, per2, CRY1, CRY2, cry1, cry2,
PER1CRY1, PER1CRY2, PER2CRY1, PER2CRY2,
REVERB, reverb, RORC, rorc, CLK, clk, BMAL1, bmal1, CLKBMAL1
```

The Julia equations in `src/CircadianFVS.jl` are a direct translation of the supplied MATLAB mammalian-circadian model.

The following internal aliases are used only inside the differential equations:

| Graph/state name | Equation alias |
|---|---|
| `PER1CRY1` | `PC11` |
| `PER1CRY2` | `PC12` |
| `PER2CRY1` | `PC21` |
| `PER2CRY2` | `PC22` |
| `CLKBMAL1` | `CB1` |

The desired attractor used in the original MATLAB code is

```text
[0.1229, 0.4361, 0.1040, 0.5496, 2.1465, 1.4523, 0.1700,
 5.0701, 0.7137, 0.3011, 2.6030, 1.0467, 1.4131, 0.2738,
 0.0971, 0.2326, 0.6503, 0.3251, 1.2300, 0.7011, 1.6326]
```

Because these values were supplied to four decimal places, the residual `f(x*)` is small but not exactly zero numerically.

## Feedback vertex sets

The six FVSs carried over from the MATLAB code are:

```text
FVS 1 = PER1, PER2, CRY1, CRY2, CLK, BMAL1, RORC
FVS 2 = PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, RORC
FVS 3 = PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, RORC
FVS 4 = PER1, PER2, CRY1, CRY2, CLK, BMAL1, rorc
FVS 5 = PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, rorc
FVS 6 = PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, rorc
```

The Julia code independently checks each candidate: after deleting those vertices, it runs a topological-sort test on the remaining directed graph. A successful topological ordering confirms that the remaining graph is acyclic.

## Numerical control experiment

For a selected FVS `S`, the experiment starts from

```julia
x0 = zeros(21)
```

and replaces only the FVS coordinates by their desired attractor values:

```julia
x0[S] = x_target[S]
```

During simulation, the controlled coordinates satisfy

```text
dx_i/dt = 0,     i in S,
```

so they remain fixed at the desired attractor values. Every other state continues to obey its original nonlinear differential equation.

This is exactly the type of override used in the supplied MATLAB script.

The Julia implementation uses a fixed-step fourth-order Runge-Kutta method. The default experiment uses

```text
final time = 25
time step  = 0.001
tolerance  = 1%
```

to stay close to the sampling used in the MATLAB checking script.

## Repository structure

```text
mammalian-circadian-fvs-julia/
├── README.md
├── Project.toml
├── run_experiment.jl
├── Circadian_FVS_Colab.ipynb
├── src/
│   └── CircadianFVS.jl
└── results/
```

### `src/CircadianFVS.jl`

Contains:

- the 21-state nonlinear model;
- the desired attractor;
- the six FVSs;
- the directed influence structure;
- state-name/index conversion;
- FVS verification by acyclicity after deletion;
- FVS-controlled dynamics;
- a Runge-Kutta simulator;
- convergence and relative-error calculations.

### `run_experiment.jl`

Runs all six FVS experiments and produces:

```text
results/summary.csv
results/max_relative_error.png
results/fvs2_normalized_states.png
```

`max_relative_error.png` compares the maximum relative state error for all six FVS choices.

`fvs2_normalized_states.png` shows the uncontrolled states for FVS 2 after normalizing each state by its target value. Convergence to `1` means convergence to the desired attractor.

### `Circadian_FVS_Colab.ipynb`

A self-contained Google Colab notebook. It installs Julia and `Plots.jl`, writes the same Julia source files into the Colab runtime, executes the complete experiment, and displays the generated plots.

## Run in Google Colab

Upload `Circadian_FVS_Colab.ipynb` to Google Colab and run the cells from top to bottom.

No MATLAB installation is required.

The notebook uses a standard Colab Python runtime only to bootstrap Julia. The dynamical model, FVS verification, numerical integration, convergence calculation, and plot generation are executed by Julia.

## Run locally with Julia

Install Julia 1.10 or later, enter the repository directory, and run

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run_experiment.jl
```

## Changing the controlled FVS

In Julia:

```julia
include("src/CircadianFVS.jl")
using .CircadianFVS

fvs = FVS_SETS[1]
t, X = rk4_simulate(fvs; x0=zeros(21), tf=25.0, dt=1e-3)
```

To verify a candidate FVS:

```julia
is_fvs, topological_order = verify_fvs(fvs)
```

## Extending from a fixed attractor to a time-varying desired trajectory

The original MATLAB implementation clamps the selected FVS states to constant attractor values. The same architecture can be generalized to a time-varying desired natural trajectory `x_target(t)` by replacing the constant clamp with

```text
x_i(t) = x_target_i(t),     i in FVS.
```

The non-FVS equations can then be integrated while reading the controlled coordinates from that reference trajectory.

That extension is not enabled by default here because the uploaded experiment supplies a fixed attractor rather than a stored periodic or otherwise time-varying target trajectory.

## MATLAB-to-Julia correspondence

| MATLAB file/function | Julia counterpart |
|---|---|
| `mammcircfunc.m` | `rhs` |
| `mammcircoverride.m` | `controlled_rhs` + exact FVS clamp |
| `mammcirc_graph_init.m` | `INFLUENCERS` |
| `mammcirc_fvs_init.m` | `FVS_SETS` |
| `name_position.m` | `state_indices` |
| `tolconverge.m` | `convergence_time` |
| `mammcirccheck.m` | `run_experiment.jl` |

## Intended use

This repository is useful as a compact computational demonstration of structural/FVS-based control for a nonlinear biological network: only a subset of state variables is directly prescribed, while the rest of the network evolves according to its original nonlinear dynamics and is tested for convergence toward the desired natural attractor.
