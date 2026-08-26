include(joinpath(@__DIR__, "src", "CircadianFVS.jl"))
using .CircadianFVS
using Printf
using Plots

mkpath(joinpath(@__DIR__, "results"))

println("Feedback Vertex Set control of the 21-state mammalian circadian model")
println("====================================================================")
@printf("Maximum |f(x*)| at the supplied attractor: %.6e\n\n", maximum(abs.(rhs(ATTRACTOR))))

println("Checking the six supplied FVS candidates:")
for (i, fvs) in enumerate(FVS_SETS)
    ok, topo = verify_fvs(fvs)
    println("  FVS $i: ", ok ? "valid (remaining graph is acyclic)" : "NOT an FVS")
end
println()

tf = 25.0
dt = 1e-3
tol = 0.01
x0 = zeros(21)

all_t = Vector{Vector{Float64}}()
all_X = Vector{Matrix{Float64}}()
all_err = Vector{Vector{Float64}}()

summary_path = joinpath(@__DIR__, "results", "summary.csv")
open(summary_path, "w") do io
    println(io, "fvs,convergence_time,max_relative_error_at_final_time")
    for (i, fvs) in enumerate(FVS_SETS)
        t, X = rk4_simulate(fvs; x0=x0, tf=tf, dt=dt)
        err = max_relative_error(X)
        tc = convergence_time(t, X; tol=tol)

        push!(all_t, t)
        push!(all_X, X)
        push!(all_err, err)

        tc_text = tc === nothing ? "not reached" : @sprintf("%.4f", tc)
        @printf("FVS %d: convergence time = %-11s final max relative error = %.6f %%\n",
                i, tc_text, 100 * err[end])

        csv_tc = tc === nothing ? "" : string(tc)
        println(io, "$(i),$(csv_tc),$(err[end])")
    end
end

# Figure 1: maximum relative error for all six FVS choices.
p1 = plot(
    xlabel="Time",
    ylabel="Maximum relative error",
    title="FVS control drives the full state toward the supplied attractor",
    yscale=:log10,
    legend=:topright,
    grid=true
)
for i in eachindex(FVS_SETS)
    stride = 10
    plot!(p1, all_t[i][1:stride:end], all_err[i][1:stride:end], label="FVS $i", linewidth=2)
end
hline!(p1, [tol], label="1% tolerance", linestyle=:dash, linewidth=2)
savefig(p1, joinpath(@__DIR__, "results", "max_relative_error.png"))

# Figure 2: normalized non-FVS states for one representative set.
chosen = 2
chosen_fvs = FVS_SETS[chosen]
chosen_idx = Set(state_indices(chosen_fvs))
t = all_t[chosen]
X = all_X[chosen]

p2 = plot(
    xlabel="Time",
    ylabel="State / target value",
    title="Non-FVS states under FVS $chosen control",
    legend=:outerright,
    grid=true
)
stride = 10
for j in 1:length(STATE_NAMES)
    if !(j in chosen_idx)
        plot!(
            p2,
            t[1:stride:end],
            X[1:stride:end, j] ./ ATTRACTOR[j],
            label=STATE_NAMES[j],
            linewidth=1.8
        )
    end
end
hline!(p2, [1.0], label="desired attractor", linestyle=:dash, linewidth=2)
savefig(p2, joinpath(@__DIR__, "results", "fvs2_normalized_states.png"))

println()
println("Saved:")
println("  results/summary.csv")
println("  results/max_relative_error.png")
println("  results/fvs2_normalized_states.png")
