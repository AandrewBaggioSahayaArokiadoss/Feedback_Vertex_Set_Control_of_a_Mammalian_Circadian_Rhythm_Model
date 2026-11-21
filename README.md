# Feedback Vertex Set Control of a 21-State Mammalian Circadian Rhythm Model

A reproducible **Python/Jupyter notebook** implementing a published 21-state nonlinear mammalian circadian clock model, its directed state-dependency graph, and six idealized feedback-vertex-set (FVS) control experiments.

The notebook includes the complete 132-parameter model, an equilibrium solver, directed graph validation, numerical integration, an uncontrolled baseline, convergence diagnostics, and plots. It has been **run end-to-end** on the supplied model.

## Quick start

1. Install Python 3.10 or newer.
2. Download `circadian_fvs_notebook.ipynb` and `requirements.txt` to the same folder.
3. In a terminal, run:

   ```bash
   python -m pip install -r requirements.txt
   jupyter lab
   ```

4. Open `circadian_fvs_notebook.ipynb` and select **Run All**. Alternatively, upload the notebook to Google Colab (Python runtime), and install any missing dependencies there.
5. The notebook creates or updates the `outputs/` folder with its CSV and PNG figures.

The notebook does **not** need the older MATLAB code, Julia code, or the original hand-drawn graph to run.

## Model and scientific sources

The mammalian clock model describes 21 molecular states: proteins, corresponding mRNAs, and protein complexes. Transcription includes multiple activating and inhibitory Hill-function factors. The chosen 132 parameters are from the FVS study, **not a newly fitted biological model**.

- **Mirsky et al. (2009),** *A model of the cell-autonomous mammalian circadian clock*, PNAS. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2699375/). Original source of the biological clock model.
- **Fiedler, Mochizuki, Kurosawa, and Saito (2013),** *Dynamics and control at feedback vertex sets. I: Informative and determining nodes in regulatory networks*, Section 7.3 and Appendix 9. [Paper with equations and parameters](https://dynamics.mi.fu-berlin.de/preprints/FiedlerMochizukiKurosawaSaito%202013-03-30%20submit%20final.pdf). Source for the precise 21 equations, the 132 parameter values, and the FVS setting studied here.

**One transcription correction:** the second activating Hill term in the `rorc` mRNA equation must use **`RORC`**, not `CLKBMAL1`. This is clear from Eq. (9.8) of Fiedler et al. The original uploaded Julia code accidentally used `CLKBMAL1` twice in the activation terms. The notebook **corrects this equation and recomputes the target equilibrium**. Using the old four-decimal target with the corrected equations produces an appreciable equilibrium residual.

## What the notebook does

The notebook:

1. Defines the 21 state variables in a fixed, documented order and loads all 132 published parameters.
2. Implements the 21-dimensional ODE right-hand side with explicit activating/repressing Hill functions and protein-complex reaction terms.
3. Starts from the approximate stationary state in the older code and computes a refined positive equilibrium using `scipy.optimize.root`. It checks `max(abs(f(x_star)))` and numerically estimates the local Jacobian eigenvalues.
4. Builds the directed influence graph from a documented dependency mapping. An edge `u -> v` means that `u` appears in the equation for `dv/dt`; **self-dependencies are omitted by convention**. The graph has **21 vertices and 78 directed edges**. An independent finite-difference sensitivity check compares the dependency mapping with the equations at several positive states.
5. Tests six seven-state FVS candidates by removing their vertices and checking that the remaining graph is acyclic. It also tests inclusion-minimality. It **does not** claim that cardinality seven has been proven globally optimal.
6. Simulates **ideal constant-value clamping** of each FVS to the recomputed equilibrium while the remaining states evolve. It also simulates an uncontrolled baseline starting from the all-zero state.
7. Uses adaptive SciPy `DOP853` integration, checks numerical nonnegativity of the recorded trajectory, and reports the final maximum relative error and Euclidean distance from the equilibrium.
8. Measures when each trajectory last enters the **1% componentwise-relative-error band** over the observed time window. This is a *finite-time, sample-grid statistic*, not a mathematical proof of convergence.
9. Generates a corrected digraph, an adjacency matrix, all 21 FVS 2 state traces, distance comparisons, and a CSV summary. A separate classical RK4 function cross-checks the older Julia approach, including step-size refinement.

## Influence Digraph

The figure below shows the **influence digraph** of the 21-state mammalian circadian rhythm model.

Each vertex represents a state variable, and a directed edge from \(u\) to \(v\) indicates that the variable \(u\) appears in the differential equation governing \(v\). Self-dependencies are omitted.

The vertices are colour-coded according to their biological roles:

- **Blue:** Protein species
- **Green:** mRNA species
- **Purple:** Protein complexes

![Influence digraph of the 21-state mammalian circadian rhythm model](outputs/circadian_dependency_digraph.png)

The digraph was generated independently from the differential equations of the published model, rather than reproducing the original figure. It is used to verify the six candidate feedback vertex sets considered in the numerical simulations; the notebook does not automatically discover those sets.

The underlying mathematical model and its regulatory interactions are attributed to Mirsky et al. (2009) and the subsequent work of Fiedler et al. (2013).

### The six FVS candidates

| Set | Controlled state variables |
|---|---|
| FVS 1 | PER1, PER2, CRY1, CRY2, CLK, BMAL1, RORC |
| FVS 2 | PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, RORC |
| FVS 3 | PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, RORC |
| FVS 4 | PER1, PER2, CRY1, CRY2, CLK, BMAL1, rorc |
| FVS 5 | PER1, PER2, CRY1, CRY2, BMAL1, CLKBMAL1, rorc |
| FVS 6 | PER1, PER2, CRY1, CRY2, CLK, CLKBMAL1, rorc |

All six sets are FVSs and inclusion-minimal for the documented influence graph. They are carried forward from the earlier project code; their original derivation is not documented in this repository.

## Outputs

Running the notebook creates:

```text
outputs/
├── circadian_dependency_digraph.png
├── circadian_dependency_matrix.png
├── circadian_fvs2_21_trajectories.png
├── circadian_fvs_distance_comparison.png
├── circadian_fvs_relative_errors.png
└── circadian_fvs_summary.csv
```

The figures are generated from the model's actual dependency dictionary rather than from the old diagram. In particular, `CLKBMAL1` is represented **once**, and the corrected `RORC -> rorc` influence is present.

### Selected numerical cross-checks

On the supplied parameter set, with a simulation interval of `0` to `25`, 2,501 saved time samples, and adaptive solver tolerances `rtol=1e-9`, `atol=1e-11`:

| Numerical diagnostic | Result |
|---|---:|
| Equilibrium residual, max absolute derivative | ~`1.8e-15` |
| Uncontrolled equilibrium Jacobian, max real eigenvalue | ~`+0.0028234` |
| Final max relative error, uncontrolled | ~`115.61%` |
| Final max relative error, FVS 2 | ~`0.520%` |
| Last sampled entry into the 1% band, FVS 2 | ~`21.83` |

These are **numerical results for one specified setup**, not general guarantees. Results can vary slightly with numerical tolerances or software versions. The positive real part of an uncontrolled Jacobian eigenvalue indicates **local instability** of the target equilibrium for the uncontrolled dynamics. It is therefore more accurate to call it a **target equilibrium** than an uncontrolled attractor.

### Why RK4 and `DOP853` can disagree near zero

Some Hill exponents are fractional and below one. The zero initial condition lies on the boundary of the nonnegative state space, so the relevant nonlinear functions need not be differentiable there. Consequently the original fixed-step RK4 scheme may converge more slowly with decreasing step size than its formal fourth-order accuracy would suggest. The notebook provides a refinement check over a short interval. It uses the positive part `max(z, 0)` only for evaluating Hill terms at temporary negative adaptive-solver *trial stages*, while checking that the accepted, recorded state trajectory remains nonnegative.

## Implementation history and attribution

**Inherited:** The biological model, all 21 equations, the parameter choices, and the general FVS-control framework come from the cited publications. The six FVS candidates were available in the earlier project code.

**Original implementation:** The researcher independently implemented and simulated this published model in **MATLAB during 2021–2023**. An important numerical difficulty was identifying a meaningful stationary point among difficult-to-interpret solver outputs; supplying a useful initial guess to MATLAB's numerical/symbolic solver enabled recovery of the desired solution.

**Later conversion:** The MATLAB implementation was subsequently converted into Julia with **LLM assistance**. This notebook is a further **AI-assisted Python/Jupyter adaptation** of that Julia version. The scientific model and parameters are not presented as new discoveries, and the present Python translation should not be misrepresented as entirely handwritten without assistance.

## Limitations

- **Ideal clamping:** Controlled variables are held at the target from time zero. This is not a physical actuator or a realizable feedback controller design.
- **Finite-time simulation:** A trajectory entering a tolerance band and staying there until time 25 is not proof of global or asymptotic convergence.
- **Equilibrium, not a globally attracting state:** The selected target is an unstable stationary state of the *uncontrolled* dynamics for these parameters. A numerical nonlinear solver only finds a stationary solution near its initial guess.
- **Graph convention:** The influence graph records off-diagonal dependencies, not the sign, magnitude, or time-varying strength of interactions. An edge does not necessarily mean activation.
- **Parameter and reproduction scope:** Use the **Fiedler et al.** parameter values if comparing results to this notebook; the original Mirsky et al. paper also discusses other parameterizations and dynamics.

## Repository layout

```text
circadian_fvs_notebook.ipynb   # Executable and annotated Python notebook
README.md                      # This document
requirements.txt               # Python dependencies
outputs/                       # Generated figures and results; commit at least circadian_dependency_digraph.png
```
