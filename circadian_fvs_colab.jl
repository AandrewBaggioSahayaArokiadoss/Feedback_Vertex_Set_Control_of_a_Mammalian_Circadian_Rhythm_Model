# Mammalian Circadian Rhythm Control Using Feedback Vertex Sets
# -----------------------------------------------------------------------------
# Standalone Julia script for Google Colab (Julia runtime already available).
#
# What it does:
#   1. Implements the 21-state mammalian circadian-rhythm ODE model.
#   2. Builds the state-influence digraph used in the supplied MATLAB files.
#   3. Verifies the six supplied feedback vertex sets (FVSs).
#   4. For each FVS, initializes those states at the desired attractor and
#      clamps them there (their derivatives are set to zero), exactly matching
#      the control interpretation used in the supplied MATLAB code.
#   5. Simulates the remaining states and checks convergence to the attractor.
#   6. Displays and saves summary plots and a CSV table.
#
# Run in Colab with:
#     include("circadian_fvs_colab.jl")
# -----------------------------------------------------------------------------

module CircadianFVS

export STATE_NAMES, ATTRACTOR, FVS_SETS, INFLUENCERS,
       state_indices, circadian_rhs, circadian_controlled_rhs, rk4_simulate,
       convergence_time, max_relative_error, verify_fvs

const STATE_NAMES = [
    "PER1", "PER2", "per1", "per2", "CRY1", "CRY2", "cry1", "cry2",
    "PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2",
    "REVERB", "reverb", "RORC", "rorc", "CLK", "clk", "BMAL1", "bmal1",
    "CLKBMAL1"
]

const STATE_INDEX = Dict(name => i for (i, name) in enumerate(STATE_NAMES))

const ATTRACTOR = [
    0.1229, 0.4361, 0.1040, 0.5496, 2.1465, 1.4523, 0.1700, 5.0701,
    0.7137, 0.3011, 2.6030, 1.0467, 1.4131, 0.2738, 0.0971, 0.2326,
    0.6503, 0.3251, 1.2300, 0.7011, 1.6326
]

const FVS_SETS = [
    ["PER1", "PER2", "CRY1", "CRY2", "CLK", "BMAL1", "RORC"],
    ["PER1", "PER2", "CRY1", "CRY2", "BMAL1", "CLKBMAL1", "RORC"],
    ["PER1", "PER2", "CRY1", "CRY2", "CLK", "CLKBMAL1", "RORC"],
    ["PER1", "PER2", "CRY1", "CRY2", "CLK", "BMAL1", "rorc"],
    ["PER1", "PER2", "CRY1", "CRY2", "BMAL1", "CLKBMAL1", "rorc"],
    ["PER1", "PER2", "CRY1", "CRY2", "CLK", "CLKBMAL1", "rorc"]
]

# For each target state, list the states that influence its differential equation.
# Edges in the interaction digraph are therefore influencer -> target.
const INFLUENCERS = Dict(
    "PER1"      => ["PER1CRY1", "PER1CRY2", "per1", "CRY1", "CRY2"],
    "PER2"      => ["PER2CRY1", "PER2CRY2", "per2", "CRY1", "CRY2"],
    "per1"      => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1"],
    "per2"      => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1"],
    "CRY1"      => ["PER1", "PER2", "cry1", "PER1CRY1", "PER2CRY1"],
    "CRY2"      => ["PER1", "PER2", "cry2", "PER1CRY2", "PER2CRY2"],
    "cry1"      => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1", "REVERB", "RORC"],
    "cry2"      => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1", "REVERB", "RORC"],
    "PER1CRY1"  => ["PER1", "CRY1"],
    "PER1CRY2"  => ["PER1", "CRY2"],
    "PER2CRY1"  => ["PER2", "CRY1"],
    "PER2CRY2"  => ["PER2", "CRY2"],
    "REVERB"    => ["reverb"],
    "reverb"    => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1"],
    "RORC"      => ["rorc"],
    "rorc"      => ["PER1CRY1", "PER1CRY2", "PER2CRY1", "PER2CRY2", "CLKBMAL1", "REVERB", "RORC"],
    "CLK"       => ["CLKBMAL1", "BMAL1", "clk"],
    "clk"       => ["REVERB", "RORC"],
    "BMAL1"     => ["CLKBMAL1", "CLK", "bmal1"],
    "bmal1"     => ["REVERB", "RORC"],
    "CLKBMAL1"  => ["CLK", "BMAL1"]
)

function state_indices(subset)
    idx = Int[]
    for name in subset
        haskey(STATE_INDEX, name) || error("Unknown state name: $name")
        push!(idx, STATE_INDEX[name])
    end
    idx
end

"""
    circadian_rhs(x)

Right-hand side of the 21-state mammalian circadian rhythm model translated
from the supplied MATLAB implementation.

State aliases used internally:
- PER1CRY1 -> PC11
- PER1CRY2 -> PC12
- PER2CRY1 -> PC21
- PER2CRY2 -> PC22
- CLKBMAL1 -> CB1
"""
function circadian_rhs(x::AbstractVector{<:Real})
    length(x) == 21 || error("Expected a 21-state vector.")

    (PER1, PER2, per1, per2, CRY1, CRY2, cry1, cry2,
     PC11, PC12, PC21, PC22, REVERB, reverb, RORC, rorc,
     CLK, clk, BMAL1, bmal1, CB1) = x

    # per1 parameters
    v0_per1 = 0.000001
    v1_per1 = 3.0
    KA1_per1 = 1.98
    na1_per1 = 2.0
    KI1_per1 = 1.07
    KI2_per1 = 3.96
    KI3_per1 = 1.68
    KI4_per1 = 3.11
    ni1_per1 = 2.0
    ni2_per1 = 1.0
    ni3_per1 = 2.0
    ni4_per1 = 4.0
    km_per1 = 2.18

    # per2 parameters
    v0_per2 = 0.09
    v1_per2 = 3.29
    KA1_per2 = 1.90
    na1_per2 = 10.0
    KI1_per2 = 4.51
    KI2_per2 = 2.98
    KI3_per2 = 2.24
    KI4_per2 = 3.31
    ni1_per2 = 1.0
    ni2_per2 = 1.0
    ni3_per2 = 9.0
    ni4_per2 = 8.0
    km_per2 = 0.20

    # cry1 parameters
    v0_cry1 = 0.26
    v1_cry1 = 2.44
    v2_cry1 = 2.89
    KA1_cry1 = 1.46
    KA2_cry1 = 3.76
    na1_cry1 = 4.91
    na2_cry1 = 3.01
    KI1_cry1 = 0.03
    KI2_cry1 = 0.77
    KI3_cry1 = 3.59
    KI4_cry1 = 3.44
    KI5_cry1 = 2.82
    ni1_cry1 = 1.0
    ni2_cry1 = 1.0
    ni3_cry1 = 6.0
    ni4_cry1 = 4.0
    ni5_cry1 = 2.24
    km_cry1 = 0.22

    # cry2 parameters
    v0_cry2 = 1.29
    v1_cry2 = 2.72
    v2_cry2 = 0.10
    KA1_cry2 = 0.69
    KA2_cry2 = 2.96
    na1_cry2 = 4.39
    na2_cry2 = 4.43
    KI1_cry2 = 4.63
    KI2_cry2 = 2.95
    KI3_cry2 = 3.57
    KI4_cry2 = 2.75
    KI5_cry2 = 3.97
    ni1_cry2 = 1.0
    ni2_cry2 = 1.0
    ni3_cry2 = 4.0
    ni4_cry2 = 8.0
    ni5_cry2 = 1.75
    km_cry2 = 0.41

    # reverb parameters
    v1_reverb = 11.06
    KA1_reverb = 3.15
    na1_reverb = 4.40
    KI1_reverb = 3.56
    KI2_reverb = 3.62
    KI3_reverb = 4.71
    KI4_reverb = 1.23
    ni1_reverb = 0.15
    ni2_reverb = 0.30
    ni3_reverb = 7.0
    ni4_reverb = 7.0
    km_reverb = 0.60

    # clk parameters
    v0_clk = 3.98
    v1_clk = 3.36
    KA1_clk = 1.59
    na1_clk = 3.50
    KI1_clk = 0.83
    ni1_clk = 1.96
    km_clk = 3.19

    # bmal1 parameters
    v0_bmal1 = 1.98
    v1_bmal1 = 4.12
    KA1_bmal1 = 2.59
    na1_bmal1 = 4.13
    KI1_bmal1 = 2.47
    ni1_bmal1 = 0.02
    km_bmal1 = 1.42

    # rorc parameters
    v0_rorc = 0.06
    v1_rorc = 3.55
    v2_rorc = 0.46
    KA1_rorc = 4.30
    KA2_rorc = 4.89
    na1_rorc = 1.57
    na2_rorc = 0.56
    KI1_rorc = 3.49
    KI2_rorc = 2.34
    KI3_rorc = 2.71
    KI4_rorc = 2.09
    KI5_rorc = 3.36
    ni1_rorc = 1.0
    ni2_rorc = 1.0
    ni3_rorc = 7.0
    ni4_rorc = 7.0
    ni5_rorc = 4.33
    km_rorc = 1.50

    # protein / complex parameters
    t_per1 = 3.05
    kp_PER1 = 2.58
    t_per2 = 2.38
    kp_PER2 = 3.0
    t_cry1 = 3.94
    kp_CRY1 = 0.312
    t_cry2 = 1.69
    kp_CRY2 = 5.9
    t_reverb = 1.60
    kp_REVERB = 0.31
    t_clk = 3.04
    kp_CLK = 1.52
    kp_BMAL1 = 2.28
    t_bmal1 = 4.00
    kp_RORC = 3.33
    t_rorc = 1.39

    # common parameters
    a_CB1 = 1.98
    d_CB1 = 0.97
    a_PC11 = 3.57
    a_PC12 = 3.12
    a_PC21 = 3.81
    a_PC22 = 4.0
    d_PC11 = 1.32
    d_PC12 = 1.85
    d_PC21 = 1.37
    d_PC22 = 2.42

    PER1_dot = t_per1 * per1 - a_PC11 * PER1 * CRY1 - a_PC12 * PER1 * CRY2 +
               d_PC11 * PC11 + d_PC12 * PC12 - kp_PER1 * PER1

    PER2_dot = t_per2 * per2 - a_PC21 * PER2 * CRY1 - a_PC22 * PER2 * CRY2 +
               d_PC21 * PC21 + d_PC22 * PC22 - kp_PER2 * PER2

    per1h1 = CB1^na1_per1 / (CB1^na1_per1 + KA1_per1^na1_per1)
    per1hi1 = KI1_per1^ni1_per1 / (KI1_per1^ni1_per1 + PC11^ni1_per1)
    per1hi2 = KI2_per1^ni2_per1 / (KI2_per1^ni2_per1 + PC12^ni2_per1)
    per1hi3 = KI3_per1^ni3_per1 / (KI3_per1^ni3_per1 + PC21^ni3_per1)
    per1hi4 = KI4_per1^ni4_per1 / (KI4_per1^ni4_per1 + PC22^ni4_per1)
    per1_dot = (v0_per1 + v1_per1 * per1h1) * per1hi1 * per1hi2 * per1hi3 * per1hi4 -
               km_per1 * per1

    per2h1 = CB1^na1_per2 / (CB1^na1_per2 + KA1_per2^na1_per2)
    per2hi1 = KI1_per2^ni1_per2 / (KI1_per2^ni1_per2 + PC11^ni1_per2)
    per2hi2 = KI2_per2^ni2_per2 / (KI2_per2^ni2_per2 + PC12^ni2_per2)
    per2hi3 = KI3_per2^ni3_per2 / (KI3_per2^ni3_per2 + PC21^ni3_per2)
    per2hi4 = KI4_per2^ni4_per2 / (KI4_per2^ni4_per2 + PC22^ni4_per2)
    per2_dot = (v0_per2 + v1_per2 * per2h1) * per2hi1 * per2hi2 * per2hi3 * per2hi4 -
               km_per2 * per2

    CRY1_dot = t_cry1 * cry1 - a_PC11 * PER1 * CRY1 - a_PC21 * PER2 * CRY1 +
               d_PC11 * PC11 + d_PC21 * PC21 - kp_CRY1 * CRY1

    CRY2_dot = t_cry2 * cry2 - a_PC12 * PER1 * CRY2 - a_PC22 * PER2 * CRY2 +
               d_PC12 * PC12 + d_PC22 * PC22 - kp_CRY2 * CRY2

    cry1h1 = CB1^na1_cry1 / (CB1^na1_cry1 + KA1_cry1^na1_cry1)
    cry1h2 = RORC^na2_cry1 / (RORC^na2_cry1 + KA2_cry1^na2_cry1)
    cry1hi1 = KI1_cry1^ni1_cry1 / (KI1_cry1^ni1_cry1 + PC11^ni1_cry1)
    cry1hi2 = KI2_cry1^ni2_cry1 / (KI2_cry1^ni2_cry1 + PC12^ni2_cry1)
    cry1hi3 = KI3_cry1^ni3_cry1 / (KI3_cry1^ni3_cry1 + PC21^ni3_cry1)
    cry1hi4 = KI4_cry1^ni4_cry1 / (KI4_cry1^ni4_cry1 + PC22^ni4_cry1)
    cry1hi5 = KI5_cry1^ni5_cry1 / (KI5_cry1^ni5_cry1 + REVERB^ni5_cry1)
    cry1_dot = (v0_cry1 + v1_cry1 * cry1h1 + v2_cry1 * cry1h2) *
               cry1hi1 * cry1hi2 * cry1hi3 * cry1hi4 * cry1hi5 - km_cry1 * cry1

    cry2h1 = CB1^na1_cry2 / (CB1^na1_cry2 + KA1_cry2^na1_cry2)
    cry2h2 = RORC^na2_cry2 / (RORC^na2_cry2 + KA2_cry2^na2_cry2)
    cry2hi1 = KI1_cry2^ni1_cry2 / (KI1_cry2^ni1_cry2 + PC11^ni1_cry2)
    cry2hi2 = KI2_cry2^ni2_cry2 / (KI2_cry2^ni2_cry2 + PC12^ni2_cry2)
    cry2hi3 = KI3_cry2^ni3_cry2 / (KI3_cry2^ni3_cry2 + PC21^ni3_cry2)
    cry2hi4 = KI4_cry2^ni4_cry2 / (KI4_cry2^ni4_cry2 + PC22^ni4_cry2)
    cry2hi5 = KI5_cry2^ni5_cry2 / (KI5_cry2^ni5_cry2 + REVERB^ni5_cry2)
    cry2_dot = (v0_cry2 + v1_cry2 * cry2h1 + v2_cry2 * cry2h2) *
               cry2hi1 * cry2hi2 * cry2hi3 * cry2hi4 * cry2hi5 - km_cry2 * cry2

    PC11_dot = a_PC11 * PER1 * CRY1 - d_PC11 * PC11
    PC12_dot = a_PC12 * PER1 * CRY2 - d_PC12 * PC12
    PC21_dot = a_PC21 * PER2 * CRY1 - d_PC21 * PC21
    PC22_dot = a_PC22 * PER2 * CRY2 - d_PC22 * PC22

    REVERB_dot = t_reverb * reverb - kp_REVERB * REVERB

    reverbh1 = CB1^na1_reverb / (CB1^na1_reverb + KA1_reverb^na1_reverb)
    reverbhi1 = KI1_reverb^ni1_reverb / (KI1_reverb^ni1_reverb + PC11^ni1_reverb)
    reverbhi2 = KI2_reverb^ni2_reverb / (KI2_reverb^ni2_reverb + PC12^ni2_reverb)
    reverbhi3 = KI3_reverb^ni3_reverb / (KI3_reverb^ni3_reverb + PC21^ni3_reverb)
    reverbhi4 = KI4_reverb^ni4_reverb / (KI4_reverb^ni4_reverb + PC22^ni4_reverb)
    reverb_dot = v1_reverb * reverbh1 * reverbhi1 * reverbhi2 * reverbhi3 * reverbhi4 -
                 km_reverb * reverb

    RORC_dot = t_rorc * rorc - kp_RORC * RORC

    rorch1 = CB1^na1_rorc / (CB1^na1_rorc + KA1_rorc^na1_rorc)
    rorch2 = CB1^na2_rorc / (CB1^na2_rorc + KA2_rorc^na2_rorc)
    rorchi1 = KI1_rorc^ni1_rorc / (KI1_rorc^ni1_rorc + PC11^ni1_rorc)
    rorchi2 = KI2_rorc^ni2_rorc / (KI2_rorc^ni2_rorc + PC12^ni2_rorc)
    rorchi3 = KI3_rorc^ni3_rorc / (KI3_rorc^ni3_rorc + PC21^ni3_rorc)
    rorchi4 = KI4_rorc^ni4_rorc / (KI4_rorc^ni4_rorc + PC22^ni4_rorc)
    rorchi5 = KI5_rorc^ni5_rorc / (KI5_rorc^ni5_rorc + REVERB^ni5_rorc)
    rorc_dot = (v0_rorc + v1_rorc * rorch1 + v2_rorc * rorch2) *
               rorchi1 * rorchi2 * rorchi3 * rorchi4 * rorchi5 - km_rorc * rorc

    CLK_dot = t_clk * clk - a_CB1 * CLK * BMAL1 + d_CB1 * CB1 - kp_CLK * CLK

    clkh1 = RORC^na1_clk / (RORC^na1_clk + KA1_clk^na1_clk)
    clkhi1 = KI1_clk^ni1_clk / (KI1_clk^ni1_clk + REVERB^ni1_clk)
    clk_dot = (v0_clk + v1_clk * clkh1) * clkhi1 - km_clk * clk

    CB1_dot = a_CB1 * CLK * BMAL1 - d_CB1 * CB1

    BMAL1_dot = t_bmal1 * bmal1 - a_CB1 * CLK * BMAL1 + d_CB1 * CB1 -
                kp_BMAL1 * BMAL1

    bmal1h1 = RORC^na1_bmal1 / (RORC^na1_bmal1 + KA1_bmal1^na1_bmal1)
    bmal1hi1 = KI1_bmal1^ni1_bmal1 / (KI1_bmal1^ni1_bmal1 + REVERB^ni1_bmal1)
    bmal1_dot = (v0_bmal1 + v1_bmal1 * bmal1h1) * bmal1hi1 - km_bmal1 * bmal1

    Float64[
        PER1_dot, PER2_dot, per1_dot, per2_dot, CRY1_dot, CRY2_dot,
        cry1_dot, cry2_dot, PC11_dot, PC12_dot, PC21_dot, PC22_dot,
        REVERB_dot, reverb_dot, RORC_dot, rorc_dot, CLK_dot, clk_dot,
        BMAL1_dot, bmal1_dot, CB1_dot
    ]
end

"""
    circadian_controlled_rhs(x, fvs_idx)

Evaluate the model and set the derivatives of controlled FVS states to zero.
This mirrors the supplied MATLAB `mammcircoverride` implementation.
"""
function circadian_controlled_rhs(x::AbstractVector{<:Real}, fvs_idx::AbstractVector{<:Integer})
    dx = circadian_rhs(x)
    dx[fvs_idx] .= 0.0
    dx
end

"""
    rk4_simulate(fvs; x0=zeros(21), attractor=ATTRACTOR, tf=25.0, dt=1e-3)

Fixed-step fourth-order Runge-Kutta simulation. FVS states are first set to
their desired attractor values and then kept exactly clamped there.
"""
function rk4_simulate(
    fvs = String[];
    x0 = zeros(21),
    attractor = ATTRACTOR,
    tf::Real = 25.0,
    dt::Real = 1e-3
)
    length(x0) == 21 || error("x0 must have 21 entries.")
    length(attractor) == 21 || error("attractor must have 21 entries.")
    tf > 0 || error("tf must be positive.")
    dt > 0 || error("dt must be positive.")

    fvs_idx = state_indices(fvs)
    x = Float64.(copy(x0))
    target = Float64.(attractor)

    if !isempty(fvs_idx)
        x[fvs_idx] .= target[fvs_idx]
    end

    nsteps = Int(round(tf / dt))
    t = collect(range(0.0; step = Float64(dt), length = nsteps + 1))
    X = Matrix{Float64}(undef, nsteps + 1, 21)
    X[1, :] .= x

    ffun = if isempty(fvs_idx)
        circadian_rhs
    else
        z -> circadian_controlled_rhs(z, fvs_idx)
    end

    for k in 1:nsteps
        k1 = ffun(x)
        k2 = ffun(x .+ 0.5 * dt .* k1)
        k3 = ffun(x .+ 0.5 * dt .* k2)
        k4 = ffun(x .+ dt .* k3)

        x .= x .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)

        # Re-impose the exact desired FVS trajectory. In the supplied MATLAB
        # experiment the desired behaviour is a fixed equilibrium, so this is
        # a constant clamp rather than a time-varying reference.
        if !isempty(fvs_idx)
            x[fvs_idx] .= target[fvs_idx]
        end

        X[k + 1, :] .= x
    end

    return t, X
end

"""
    max_relative_error(X, target=ATTRACTOR)

Maximum componentwise relative error at every time sample.
"""
function max_relative_error(X::AbstractMatrix{<:Real}, target = ATTRACTOR)
    size(X, 2) == length(target) || error("X and target dimensions do not match.")
    err = Vector{Float64}(undef, size(X, 1))
    for k in axes(X, 1)
        err[k] = maximum(abs.((X[k, :] .- target) ./ target))
    end
    err
end

"""
    convergence_time(t, X, target=ATTRACTOR; tol=0.01)

Earliest time after which all states remain within `tol` relative error of the
target for the remainder of the simulation. Returns `nothing` if that does not
happen within the simulated interval.
"""
function convergence_time(t, X, target = ATTRACTOR; tol::Real = 0.01)
    size(X, 1) == length(t) || error("t and X lengths do not match.")
    size(X, 2) == length(target) || error("X and target dimensions do not match.")

    inside = Vector{Bool}(undef, length(t))
    for k in eachindex(t)
        inside[k] = all(abs.((X[k, :] .- target) ./ target) .<= tol)
    end

    last_outside = findlast(!, inside)
    if last_outside === nothing
        return t[1]
    elseif last_outside == length(t)
        return nothing
    else
        return t[last_outside + 1]
    end
end

"""
    verify_fvs(fvs)

Remove the proposed FVS from the influence digraph and apply Kahn's algorithm.
Returns `(is_feedback_vertex_set, topological_order_of_remaining_graph)`.
"""
function verify_fvs(fvs)
    removed = Set(String.(fvs))
    remaining = [v for v in STATE_NAMES if !(v in removed)]

    indegree = Dict(v => 0 for v in remaining)
    outgoing = Dict(v => String[] for v in remaining)

    for target in remaining
        for source in INFLUENCERS[target]
            if !(source in removed)
                push!(outgoing[source], target)
                indegree[target] += 1
            end
        end
    end

    queue = sort([v for v in remaining if indegree[v] == 0])
    order = String[]

    while !isempty(queue)
        v = popfirst!(queue)
        push!(order, v)
        for w in outgoing[v]
            indegree[w] -= 1
            if indegree[w] == 0
                push!(queue, w)
                sort!(queue)
            end
        end
    end

    return length(order) == length(remaining), order
end

end # module


# =============================================================================
# COLAB / STANDALONE EXPERIMENT
# =============================================================================

using Printf

# Julia is assumed to already be installed. Install only Plots if the current
# Julia environment does not have it yet.
try
    @eval using Plots
catch
    import Pkg
    Pkg.add("Plots")
    @eval using Plots
end

using .CircadianFVS

println("\nFeedback Vertex Set control of the 21-state mammalian circadian model")
println("====================================================================")
@printf("Maximum |f(x*)| at the supplied (rounded) attractor: %.6e\n\n",
        maximum(abs.(CircadianFVS.circadian_rhs(ATTRACTOR))))

println("Checking the six supplied FVS candidates:")
for (i, fvs) in enumerate(FVS_SETS)
    ok, topo = verify_fvs(fvs)
    @printf("  FVS %d: %s; remaining states = %d\n",
            i,
            ok ? "VALID (remaining digraph is acyclic)" : "NOT an FVS",
            length(topo))
end
println()

# Simulation settings taken from the MATLAB verification script.
tf = 25.0
dt = 1e-3
tol = 0.01
x0 = zeros(21)

all_t = Vector{Vector{Float64}}()
all_X = Vector{Matrix{Float64}}()
all_err = Vector{Vector{Float64}}()
convergence_times = Union{Nothing,Float64}[]

println("Running the six FVS-control simulations ...")
for (i, fvs) in enumerate(FVS_SETS)
    t, X = rk4_simulate(fvs; x0=x0, tf=tf, dt=dt)
    err = max_relative_error(X)
    tc = convergence_time(t, X; tol=tol)

    push!(all_t, t)
    push!(all_X, X)
    push!(all_err, err)
    push!(convergence_times, tc)

    tc_text = tc === nothing ? "not reached" : @sprintf("%.4f", tc)
    @printf("  FVS %d: convergence time = %-11s | final max relative error = %.6f %%\n",
            i, tc_text, 100 * err[end])
end

# -----------------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------------
println("\nSummary")
println("-------")
println(rpad("FVS", 8), rpad("Convergence time", 22), "Final max relative error")
for i in eachindex(FVS_SETS)
    tc_text = convergence_times[i] === nothing ? "not reached" : @sprintf("%.4f", convergence_times[i])
    err_text = @sprintf("%.6f %%", 100 * all_err[i][end])
    println(rpad("FVS $i", 8), rpad(tc_text, 22), err_text)
end

# Save a compact CSV as well.
summary_csv = "circadian_fvs_summary.csv"
open(summary_csv, "w") do io
    println(io, "fvs,is_valid,convergence_time,final_max_relative_error")
    for i in eachindex(FVS_SETS)
        ok, _ = verify_fvs(FVS_SETS[i])
        tc_text = convergence_times[i] === nothing ? "" : string(convergence_times[i])
        println(io, "$(i),$(ok),$(tc_text),$(all_err[i][end])")
    end
end

# -----------------------------------------------------------------------------
# Plot 1: representative FVS (FVS 2), showing only the non-FVS states.
# Each state is normalized by its desired attractor value, so convergence is
# visually indicated by all curves approaching 1.
# -----------------------------------------------------------------------------
stride = 10
chosen = 2
chosen_fvs = FVS_SETS[chosen]
chosen_idx = Set(state_indices(chosen_fvs))
t = all_t[chosen]
X = all_X[chosen]

p2 = plot(
    xlabel = "Time",
    ylabel = "State / target value",
    title = "Non-FVS states under FVS $chosen control",
    legend = :outerright,
    grid = true,
    size = (1050, 650)
)

for j in 1:length(STATE_NAMES)
    if !(j in chosen_idx)
        plot!(p2,
              t[1:stride:end],
              X[1:stride:end, j] ./ ATTRACTOR[j],
              label = STATE_NAMES[j],
              linewidth = 1.8)
    end
end
hline!(p2, [1.0], label = "desired attractor", linestyle = :dash, linewidth = 2)
display(p2)
savefig(p2, "circadian_fvs2_normalized_nonfvs_states.png")

# -----------------------------------------------------------------------------
# Plot 2: individual state trajectories for the representative FVS.
# Controlled states are shown too; they stay exactly at their target values.
# -----------------------------------------------------------------------------
p3 = plot(
    layout = (7, 3),
    size = (1200, 1500),
    legend = false,
    plot_title = "All 21 states under FVS $chosen control"
)

for j in 1:21
    plot!(p3,
          t[1:stride:end],
          X[1:stride:end, j],
          subplot = j,
          linewidth = 1.5,
          title = STATE_NAMES[j],
          xlabel = "t")
    hline!(p3, [ATTRACTOR[j]], subplot = j, linestyle = :dash, linewidth = 1.2)
end

display(p3)
savefig(p3, "circadian_fvs2_all_state_trajectories.png")

println("\nGenerated files:")
println("  $summary_csv")
println("  circadian_fvs2_normalized_nonfvs_states.png")
println("  circadian_fvs2_all_state_trajectories.png")
println("\nFinished.")
