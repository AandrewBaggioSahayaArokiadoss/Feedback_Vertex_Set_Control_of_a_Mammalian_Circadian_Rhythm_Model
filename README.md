# Feedback Vertex Set control of the 21-state mammalian circadian model

## About
This Julia script simulates feedback-vertex-set (FVS) control of a 21-state mammalian circadian rhythm model. It verifies six candidate FVS choices, runs the controlled nonlinear dynamics, reports convergence statistics, and produces two visual outputs for a representative FVS (FVS 2):

1. **All 21 state trajectories** versus time.
2. **State error trajectories** versus time, where the error is defined as
   $e_i(t) = x_i(t) - x_i^*$,
   with $x_i^*$ equal to the desired attractor value for state $i$.

In the controlled simulations, the chosen FVS states are clamped to their target attractor values, while the remaining states evolve according to the model dynamics.

## What was changed
Compared with the earlier version of the script:
- the **normalized non-FVS-only plot** was removed;
- the **all-21-states plot** was retained;
- a new **21-panel error plot** was added, showing $x_i(t) - x_i^*$ for every state.

## Files generated
Running the script creates:
- `circadian_fvs_summary.csv` — convergence summary for all six FVS choices;
- `circadian_fvs2_all_state_trajectories.png` — 21 subplots of the state trajectories for FVS 2;
- `circadian_fvs2_state_errors.png` — 21 subplots of the state errors for FVS 2.

## How to run
In a Julia-enabled Colab session or any Julia environment:

```julia
include("circadian_fvs_21states_with_error.jl")
```

The script automatically installs `Plots` if it is not already available.

## Notes
- The desired long-term behavior used here is a fixed attractor, so the reference trajectory is constant.
- The error plot therefore measures deviation from that attractor, not error relative to a time-varying reference.
- The script uses **FVS 2** as the representative case for the final plots, just as in the earlier workflow.
