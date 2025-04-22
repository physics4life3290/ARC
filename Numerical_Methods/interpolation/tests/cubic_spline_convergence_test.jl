

using Plots
include("../methods/cubic_spline.jl")


function run_interp_convergence_test()
    # 1. Function the spline will ace (cubic polynomial)
    f_ace(x) = x^3 + 2x^2 + x + 1

    # 2. Function the spline will do okay on (smooth, non-polynomial)
    f_okay(x) = sin(x)

    # 3. Function that pushes the limits (discontinuous step function)
    f_limit(x) = x < 0.5 ? 0.0 : 1.0

    lengths = [10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 10000]


    ind_var_exact = collect(range(0.0, 3.0, length=100))
    dep_var_exact_ace = f_ace.(ind_var_exact)
    dep_var_exact_okay = f_okay.(ind_var_exact)
    dep_var_exact_limit = f_limit.(ind_var_exact)

    # Convergence Test #

    ys = [dep_var_exact_ace, dep_var_exact_okay, dep_var_exact_limit]

    interp_ys_ace = []
    interp_ys_okay = []
    interp_ys_limit = []
    ind_var_interps = []

    for len in lengths 
        global ind_var_interp = collect(range(0.1, 2.9, length=len))
        push!(ind_var_interps, ind_var_interp)  # Store the interpolated x values
        for y in ys
            if y == dep_var_exact_ace
                cs = cubic_spline(ind_var_exact, y)  # Create the cubic spline
                push!(interp_ys_ace, cs.(ind_var_interp))  # Interpolate the values
            elseif y == dep_var_exact_okay
                cs = cubic_spline(ind_var_exact, y)  # Create the cubic spline
                push!(interp_ys_okay, cs.(ind_var_interp))  # Interpolate the values
            elseif y == dep_var_exact_limit
                cs = cubic_spline(ind_var_exact, y)  # Create the cubic spline
                push!(interp_ys_limit, cs.(ind_var_interp))  # Interpolate the values
            end
        end
    end

    
    plot(title="Cubic Spline Interpolation\nConvergence Plot", xlabel="x", ylabel="y", size=(800, 600), grid=true, xlims=(0,1), ylims=(0.0, 0.05))
    plot!(ind_var_exact, dep_var_exact_ace ./ maximum(dep_var_exact_ace), label="Exact (Cubic)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
    plot!(ind_var_interps[1], interp_ys_ace[1] ./ maximum(interp_ys_ace[1]), label="10 points", linestyle=:dash)
    plot!(ind_var_interps[2], interp_ys_ace[2] ./ maximum(interp_ys_ace[2]), label="20 points", linestyle=:dash)
    plot!(ind_var_interps[3], interp_ys_ace[3] ./ maximum(interp_ys_ace[3]), label="50 points", linestyle=:dash)
    plot!(ind_var_interps[4], interp_ys_ace[4] ./ maximum(interp_ys_ace[4]), label="100 points", linestyle=:dash)
    plot!(ind_var_interps[5], interp_ys_ace[5] ./ maximum(interp_ys_ace[5]), label="250 points", linestyle=:dash)
    plot!(ind_var_interps[6], interp_ys_ace[6] ./ maximum(interp_ys_ace[6]), label="500 points", linestyle=:dash)
    plot!(ind_var_interps[7], interp_ys_ace[7] ./ maximum(interp_ys_ace[7]), label="1000 points", linestyle=:dash)
    plot!(ind_var_interps[8], interp_ys_ace[8] ./ maximum(interp_ys_ace[8]), label="2500 points", linestyle=:dash)
    plot!(ind_var_interps[9], interp_ys_ace[9] ./ maximum(interp_ys_ace[9]), label="5000 points", linestyle=:dash)
    plot!(ind_var_interps[10], interp_ys_ace[10] ./ maximum(interp_ys_ace[10]), label="10000 points", linestyle=:dash)
    savefig("Numerical_Methods/interpolation/tests/graphs/cubic_spline_ace.pdf")

    plot(title="Cubic Spline Interpolation\nConvergence Plot", xlabel="x", ylabel="y", size=(800, 600), grid=true, xlims=(1, 2), ylims=(0.75, 1.0))
    plot!(ind_var_exact, dep_var_exact_okay ./ maximum(dep_var_exact_okay), label="Exact (Smooth)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
    plot!(ind_var_interps[1], interp_ys_okay[1] ./ maximum(interp_ys_okay[1]), label="10 points", linestyle=:dash)
    plot!(ind_var_interps[2], interp_ys_okay[2] ./ maximum(interp_ys_okay[2]), label="20 points", linestyle=:dash)
    plot!(ind_var_interps[3], interp_ys_okay[3] ./ maximum(interp_ys_okay[3]), label="50 points", linestyle=:dash)
    plot!(ind_var_interps[4], interp_ys_okay[4] ./ maximum(interp_ys_okay[4]), label="100 points", linestyle=:dash)
    plot!(ind_var_interps[5], interp_ys_okay[5] ./ maximum(interp_ys_okay[5]), label="250 points", linestyle=:dash)
    plot!(ind_var_interps[6], interp_ys_okay[6] ./ maximum(interp_ys_okay[6]), label="500 points", linestyle=:dash)
    plot!(ind_var_interps[7], interp_ys_okay[7] ./ maximum(interp_ys_okay[7]), label="1000 points", linestyle=:dash)
    plot!(ind_var_interps[8], interp_ys_okay[8] ./ maximum(interp_ys_okay[8]), label="2500 points", linestyle=:dash)
    plot!(ind_var_interps[9], interp_ys_okay[9] ./ maximum(interp_ys_okay[9]), label="5000 points", linestyle=:dash)
    plot!(ind_var_interps[10], interp_ys_okay[10] ./ maximum(interp_ys_okay[10]), label="10000 points", linestyle=:dash)
    savefig("Numerical_Methods/interpolation/tests/graphs/cubic_spline_okay.pdf")

    plot(title="Cubic Spline Interpolation\nConvergence Plot", xlabel="x", ylabel="y", size=(800, 600), grid=true, xlims=(0.4, 0.8), ylims=(0.75, 1.25))
    plot!(ind_var_exact, dep_var_exact_limit ./ maximum(dep_var_exact_limit), label="Exact (Step)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
    plot!(ind_var_interps[1], interp_ys_limit[1], label="10 points", linestyle=:dash)
    plot!(ind_var_interps[2], interp_ys_limit[2], label="20 points", linestyle=:dash)
    plot!(ind_var_interps[3], interp_ys_limit[3], label="50 points", linestyle=:dash)
    plot!(ind_var_interps[4], interp_ys_limit[4], label="100 points", linestyle=:dash)
    plot!(ind_var_interps[5], interp_ys_limit[5], label="250 points", linestyle=:dash)
    plot!(ind_var_interps[6], interp_ys_limit[6], label="500 points", linestyle=:dash)
    plot!(ind_var_interps[7], interp_ys_limit[7], label="1000 points", linestyle=:dash)
    plot!(ind_var_interps[8], interp_ys_limit[8], label="2500 points", linestyle=:dash)
    plot!(ind_var_interps[9], interp_ys_limit[9], label="5000 points", linestyle=:dash)
    plot!(ind_var_interps[10], interp_ys_limit[10], label="10000 points", linestyle=:dash)
    savefig("Numerical_Methods/interpolation/tests/graphs/cubic_spline_limit.pdf")
end

