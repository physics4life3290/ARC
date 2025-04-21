



using Plots
include("../methods/cubic_spline.jl")

function run_interp_ptp_vs_all_points_test()
    # 1. Polynomial function (should ace interpolation)
    f_ace(x) = x^3 + 2x^2 + x + 1

    # 2. Smooth sine function
    f_okay(x) = sin(x)

    # 3. Discontinuous step function
    f_limit(x) = x < 0.5 ? 0.0 : 1.0

    lengths = [10, 100, 1000, 10000]

        # Generate data
        ind_var_exact = collect(range(0.0, 3.0, length=lengths[4]))
        dep_var_exact_ace = f_ace.(ind_var_exact)
        dep_var_exact_okay = f_okay.(ind_var_exact)
        dep_var_exact_limit = f_limit.(ind_var_exact)

    for lens in lengths

        # point to point interpolation
        interp_xs = collect(range(0.1, 2.9, length=lens))  # separate resolution
        ptp_interped_ys_ace = similar(interp_xs)
        ptp_interped_ys_okay = similar(interp_xs)
        ptp_interped_ys_limit = similar(interp_xs)

        for (j, xq) in enumerate(interp_xs)
            # Find the closest index to xq
            idx = findfirst(x -> x > xq, ind_var_exact)
            if idx === nothing || idx < 2 || idx > length(ind_var_exact) - 1
                ptp_interped_ys_ace[j] = NaN  # or use linear extrapolation, etc.
                ptp_interped_ys_okay[j] = NaN
                ptp_interped_ys_limit[j] = NaN 
                continue
            end
            xs = ind_var_exact[idx-1:idx+1]
            ys = dep_var_exact_ace[idx-1:idx+1]
            global cs = cubic_spline(xs, ys)
            ptp_interped_ys_ace[j] = cs(xq)

            ys = dep_var_exact_okay[idx-1:idx+1]
            global cs = cubic_spline(xs, ys)
            ptp_interped_ys_okay[j] = cs(xq)

            ys = dep_var_exact_limit[idx-1:idx+1]
            global cs = cubic_spline(xs, ys)
            ptp_interped_ys_limit[j] = cs(xq)
        end


        cs_all_ace = cubic_spline(ind_var_exact, dep_var_exact_ace)  # Create the cubic spline
        cs_all_okay = cubic_spline(ind_var_exact, dep_var_exact_okay)  # Create the cubic spline
        cs_all_limit = cubic_spline(ind_var_exact, dep_var_exact_limit)  # Create the cubic spline
        all_interped_ys_ace = cs_all.(interp_xs)  # Interpolate the values
        all_interped_ys_okay = cs_all_okay.(interp_xs)  # Interpolate the values
        all_interped_ys_limit = cs_all_limit.(interp_xs)  # Interpolate the values

        # Plot interpolated curve
        plot(ind_var_exact, dep_var_exact_ace, label="Cubic Spline", lw=2,
            title="Cubic Spline Interpolation of \n f(x) = x^3 + 2x^2 + x + 1",
            xlabel="x", ylabel="f(x)", legend=:topright)
        plot!(interp_xs, all_interped_ys_ace, label="All Points Interpolation", lw=2, ls=:dash)
        plot!(interp_xs, ptp_interped_ys_ace, label="Point to Point Interpolation", lw=2, ls=:dot)
        savefig("Numerical_Methods/interpolation/tests/graphs/ptp_vs_all/_$(lens)_points/cubic_spline_ace_$(lens).pdf")

        plot(ind_var_exact, dep_var_exact_okay, label="Cubic Spline", lw=2,
            title="Cubic Spline Interpolation of \n f(x) = sin(x)",
            xlabel="x", ylabel="f(x)", legend=:topright)
        plot!(interp_xs, all_interped_ys_okay, label="All Points Interpolation", lw=2, ls=:dash)
        plot!(interp_xs, ptp_interped_ys_okay, label="Point to Point Interpolation", lw=2, ls=:dot)
        savefig("Numerical_Methods/interpolation/tests/graphs/ptp_vs_all/_$(lens)_points/cubic_spline_okay_$(lens).pdf")

        plot(ind_var_exact, dep_var_exact_limit, label="Cubic Spline", lw=2,
            title="Cubic Spline Interpolation of \n f(x) = step function",
            xlabel="x", ylabel="f(x)", legend=:topright, xlims=(0.4, 0.8), ylims=(0.75, 1.25))
        plot!(interp_xs, all_interped_ys_limit, label="All Points Interpolation", lw=2, ls=:dash)
        plot!(interp_xs, ptp_interped_ys_limit, label="Point to Point Interpolation", lw=2, ls=:dot)
        savefig("Numerical_Methods/interpolation/tests/graphs/ptp_vs_all/_$(lens)_points/cubic_spline_limit_$(lens).pdf")
    end
end
