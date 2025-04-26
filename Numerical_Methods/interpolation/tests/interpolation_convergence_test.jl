using Plots
include("../interpolation_include.jl")

# Define test functions and their derivatives
f_smooth(var) = var^3
df_smooth(var) = 3 * var^2

f_oscillatory(var) = sin(var)
df_oscillatory(var) = cos(var)

f_square_wave(var) = sign(sin(var))
df_square_wave(var) = cos(var) * (var > 0 ? 1 : -1)

f_exponential(var) = exp(var)
df_exponential(var) = exp(var)

f_gaussian(var) = exp(-var^2)
df_gaussian(var) = -2 * var * exp(-var^2)

f_abs(var) = abs(var)
df_abs(var) = var > 0 ? 1 : (var < 0 ? -1 : 0)

f_highfreq(var) = sin(10 * var)
df_highfreq(var) = 10 * cos(10 * var)

f_piecewise(var) = var < 5 ? var : 2 * var - 5
df_piecewise(var) = var < 5 ? 1 : 2

# Define independent variable
independent_var = collect(range(0.0, stop=10.0, length=1000))

# Compute dependent variables
smooth_dep_var = f_smooth.(independent_var)
dsmooth_dep_var = df_smooth.(independent_var)

oscillatory_dep_var = f_oscillatory.(independent_var)
doscillatory_dep_var = df_oscillatory.(independent_var)

square_wave_dep_var = f_square_wave.(independent_var)
dsquare_wave_dep_var = df_square_wave.(independent_var)

exponential_dep_var = f_exponential.(independent_var)
dexponential_dep_var = df_exponential.(independent_var)

gaussian_dep_var = f_gaussian.(independent_var)
dgaussian_dep_var = df_gaussian.(independent_var)

abs_dep_var = f_abs.(independent_var)
dabs_dep_var = df_abs.(independent_var)

highfreq_dep_var = f_highfreq.(independent_var)
dhighfreq_dep_var = df_highfreq.(independent_var)

piecewise_dep_var = f_piecewise.(independent_var)
dpiecewise_dep_var = df_piecewise.(independent_var)

# Store variable data and labels
test_variables = Dict(
    "smooth" => (smooth_dep_var, dsmooth_dep_var),
    "oscillatory" => (oscillatory_dep_var, doscillatory_dep_var),
    "square_wave" => (square_wave_dep_var, dsquare_wave_dep_var),
    "exponential" => (exponential_dep_var, dexponential_dep_var),
    "gaussian" => (gaussian_dep_var, dgaussian_dep_var),
    "abs" => (abs_dep_var, dabs_dep_var),
    "highfreq" => (highfreq_dep_var, dhighfreq_dep_var),
    "piecewise" => (piecewise_dep_var, dpiecewise_dep_var),
)


lengths = [25, 50, 75, 100, 500, 1000]

function interpolation_convergence_test(interp_func::Function; use_deriv::Bool=false)
    all_results = Dict{String, Dict{Int, Vector{Float64}}}()
    interp_xs = Dict{Int, Vector{Float64}}()

    for (label, (dep_var, ddep_var)) in test_variables
        all_results[label] = Dict{Int, Vector{Float64}}()

        for len in lengths
            interp_ind_var = collect(range(0.1, stop=9.9, length=len))  # <-- Cover full interval [0,10]            
            interp_dep_var = similar(interp_ind_var)

            for i in 1:len
                if i == 1
                    idxs = [1, 2, 3]
                elseif i == len
                    idxs = [len-2, len-1, len]
                else
                    idxs = [i-1, i, i+1]
                end

                ind_var_interval = independent_var[idxs]
                dep_var_interval = dep_var[idxs]
                ddep_var_interval = ddep_var[idxs]
                x_eval = interp_ind_var[i]

                interp_dep_var[i] = use_deriv ?
                    interp_func(ind_var_interval, dep_var_interval, ddep_var_interval, x_eval) :
                    interp_func(ind_var_interval, dep_var_interval, x_eval)
            end

            all_results[label][len] = interp_dep_var
            interp_xs[len] = interp_ind_var
        end
    end

    # Plotting
    for (label, result_dict) in all_results
        plt = plot(title="Normalized Interpolations for $label function", xlabel="x", ylabel="f(x)", legend=:bottomright)#, xlims=(0.0, 10.0), ylims=(-1.5, 1.5))

        # Plot the exact function
        plot!(plt, independent_var, test_variables[label][1]/maximum(abs.(test_variables[label][1])), label="Exact", lw=3, ls=:dash, color=:black)

        # Plot interpolations
        for len in lengths
            plot!(plt, interp_xs[len], result_dict[len]/maximum(abs.(result_dict[len])), label="n=$len", lw=2)
        end

        #display(plt)
        savefig(plt, "Numerical_Methods/interpolation/tests/graphs/convergence/interpolation_$(label).pdf")
    end
end

# Example usage with cubic_spline
interpolation_convergence_test(cubic_spline)
#interpolation_convergence_test(newton_interpolation)
#interpolation_convergence_test(hermite_interpolation; use_deriv=true)
#interpolation_convergence_test(weno5_interpolate_at)
