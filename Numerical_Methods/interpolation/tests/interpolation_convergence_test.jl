

using Plots


function run_interp_convergence_test(interp_fn::Function)
    # Functions to test
    f_ace(x) = x^3 + 2x^2 + x + 1
    f_okay(x) = sin(x)
    f_limit(x) = x < 0.5 ? 0.0 : 1.0

    lengths = [10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 10000]

    # Exact values
    ind_var_exact = collect(range(0.0, 3.0, length=100))
    dep_var_exact_ace = f_ace.(ind_var_exact)
    dep_var_exact_okay = f_okay.(ind_var_exact)
    dep_var_exact_limit = f_limit.(ind_var_exact)

    # Output storage
    interp_ys_ace = []
    interp_ys_okay = []
    interp_ys_limit = []
    ind_var_interps = []

    for len in lengths 
        ind_var_interp = collect(range(0.1, 2.9, length=len))
        push!(ind_var_interps, ind_var_interp)

        # Evaluate functions at exact x values
        y_ace = f_ace.(ind_var_exact)
        y_okay = f_okay.(ind_var_exact)
        y_limit = f_limit.(ind_var_exact)

        # Interpolate
        push!(interp_ys_ace, interp_fn(ind_var_exact, y_ace, ind_var_interp))
        push!(interp_ys_okay, interp_fn(ind_var_exact, y_okay, ind_var_interp))
        push!(interp_ys_limit, interp_fn(ind_var_exact, y_limit, ind_var_interp))
    end

    # Plotting (reused logic)
    function plot_convergence(title_str, exact, interps, fname)
        plot(title=title_str, xlabel="x", ylabel="y", size=(800, 600), grid=true)
        plot!(ind_var_exact, exact, label="Exact")
        for (i, y) in enumerate(interps)
            plot!(ind_var_interps[i], y, label="$(lengths[i]) points", linestyle=:dash)
        end
        savefig(fname)
    end

    plot_convergence("Interpolation Convergence\n(Cubic Poly)", ind_var_interps, interp_ys_ace, "Numerical_Methods/interpolation/tests/graphs/interp_ace.pdf")
    plot_convergence("Interpolation Convergence\n(Smooth)", ind_var_interps, interp_ys_okay, "Numerical_Methods/interpolation/tests/graphs/interp_okay.pdf")
    plot_convergence("Interpolation Convergence\n(Step)", ind_var_interps, interp_ys_limit, "Numerical_Methods/interpolation/tests/graphs/interp_limit.pdf")
end


