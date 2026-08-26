# Feedback Vertex Set control of the 21-state mammalian circadian model

## About
This Julia script reproduces feedback-vertex-set (FVS) control for a 21-state mammalian circadian rhythm model. It verifies six candidate feedback vertex sets, simulates the controlled nonlinear dynamics, reports convergence statistics, and visualizes how the complete system approaches the desired long-term behavior.

The selected FVS states are initialized at their desired attractor values and held there, while the remaining states evolve according to the original nonlinear model.

## Plots
The script generates two figures for the representative choice **FVS 2**.

### 1. All 21 state trajectories
`circadian_fvs2_all_state_trajectories.png`

This figure contains one subplot for each of the 21 state variables. The solid curve is the simulated state and the dashed horizontal line is its desired attractor value.

### 2. Euclidean distance from the desired trajectory
`circadian_fvs2_euclidean_distance.png`

At every time instant, the error is computed as the Euclidean distance between the full 21-dimensional system state and the desired state:

$$
d(t)=\|x(t)-x^*\|_2
=\sqrt{\sum_{i=1}^{21}\left(x_i(t)-x_i^*\right)^2}.
$$

Here, $x(t)\in\mathbb{R}^{21}$ is the simulated trajectory and $x^*$ is the desired attractor. Therefore, $d(t)=0$ means that all 21 states simultaneously coincide with the desired long-term state. A decreasing distance indicates convergence of the complete network toward the target.

Because the desired behavior used in this experiment is a fixed equilibrium, the reference trajectory is constant. The same definition extends immediately to a time-varying reference $x^*(t)$ by replacing $x^*$ with $x^*(t)$.

## Additional output
The script also creates:

- `circadian_fvs_summary.csv` — FVS validity, the original 1%-per-state convergence time, and the final Euclidean distance for all six FVS choices.

## Running in Colab
With Julia already available in the Colab runtime, upload the `.jl` file and run:

```julia
include("circadian_fvs_21states_euclidean_error.jl")
```

The plotted error is the Euclidean distance, while the reported convergence time retains the original criterion that every state must remain within 1% relative error of its target value.

The script uses `Plots.jl`; it installs `Plots` automatically only if the package is missing.
