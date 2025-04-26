using Plots
using Statistics  # For mean()
using LinearAlgebra
using Random

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

# Updated function
using LinearAlgebra
using Random
using Plots

function interpolation_convergence_test(interp_func::Function; use_deriv::Bool=false, clustered_points::Bool=false)
    all_results = Dict{String, Dict{Int, Vector{Float64}}}()
    interp_xs = Dict{Int, Vector{Float64}}()

    errors_Linf = Dict{String, Vector{Float64}}()
    errors_L2 = Dict{String, Vector{Float64}}()

    rng = MersenneTwister(1234)  # Fixed seed for reproducibility

    for (label, (dep_var, ddep_var)) in test_variables
        all_results[label] = Dict{Int, Vector{Float64}}()
        errors_Linf[label] = Float64[]
        errors_L2[label] = Float64[]

        for len in lengths
            # Random or clustered points
            if clustered_points
                interp_ind_var = sort(rand(rng, len).^3 * (9.9 - 0.1) .+ 0.1)
            else
                interp_ind_var = sort(rand(rng, len) * (9.9 - 0.1) .+ 0.1)
            end

            interp_dep_var = similar(interp_ind_var)
            true_dep_var = similar(interp_ind_var)

            for i in 1:len
                idx = searchsortedfirst(independent_var, interp_ind_var[i])

                if idx <= 2
                    idxs = [1, 2, 3]
                elseif idx >= length(independent_var) - 1
                    idxs = [length(independent_var)-2, length(independent_var)-1, length(independent_var)]
                else
                    idxs = [idx-1, idx, idx+1]
                end

                ind_var_interval = independent_var[idxs]
                dep_var_interval = dep_var[idxs]
                ddep_var_interval = ddep_var[idxs]
                x_eval = interp_ind_var[i]

                interp_dep_var[i] = use_deriv ? 
                    interp_func(ind_var_interval, dep_var_interval, ddep_var_interval, x_eval) :
                    interp_func(ind_var_interval, dep_var_interval, x_eval)

                true_dep_var[i] = f_value_from_label(label, x_eval)
            end

            all_results[label][len] = interp_dep_var
            interp_xs[len] = interp_ind_var

            # Calculate errors
            err = abs.(true_dep_var .- interp_dep_var)
            push!(errors_Linf[label], maximum(err))
            push!(errors_L2[label], sqrt(mean(err.^2)))
        end
    end

    # Plotting interpolations
    for (label, result_dict) in all_results
        plt = plot(title="Normalized Interpolations for $label function", xlabel="x", ylabel="f(x)", legend=:bottomright)
        plot!(plt, independent_var, test_variables[label][1]/maximum(abs.(test_variables[label][1])), label="Exact", lw=3, ls=:dash, color=:black)

        for len in lengths
            plot!(plt, interp_xs[len], result_dict[len]/maximum(abs.(result_dict[len])), label="n=$len", lw=2)
        end

        savefig(plt, "Numerical_Methods/interpolation/tests/graphs/interpolation_$(label).pdf")
    end

    # Plotting convergence graphs and computing rates
    for label in keys(test_variables)
        safe_Linf = max.(errors_Linf[label], 1e-15)
        safe_L2 = max.(errors_L2[label], 1e-15)

        plt_err = plot(
            xscale=:log10, 
            yscale=:log10, 
            title="Convergence plot for $label", 
            xlabel="Number of points (log)", 
            ylabel="Error (log)", 
            legend=:bottomleft
        )
        plot!(plt_err, lengths, safe_Linf, label="L∞ error", marker=:circle, lw=2)
        plot!(plt_err, lengths, safe_L2, label="L2 error", marker=:square, lw=2)
        savefig(plt_err, "Numerical_Methods/interpolation/tests/graphs/error_convergence_$(label).pdf")

        # Compute and print convergence rates
        rate_Linf = compute_convergence_rate(lengths, safe_Linf)
        rate_L2 = compute_convergence_rate(lengths, safe_L2)

        println("Convergence rates for $label: L∞ = $(round(rate_Linf, digits=2)), L2 = $(round(rate_L2, digits=2))")
    end

end

# Helper: Get true function value from label
function f_value_from_label(label::String, x::Float64)
    if label == "smooth"
        return f_smooth(x)
    elseif label == "oscillatory"
        return f_oscillatory(x)
    elseif label == "square_wave"
        return f_square_wave(x)
    elseif label == "exponential"
        return f_exponential(x)
    elseif label == "gaussian"
        return f_gaussian(x)
    elseif label == "abs"
        return f_abs(x)
    elseif label == "highfreq"
        return f_highfreq(x)
    elseif label == "piecewise"
        return f_piecewise(x)
    else
        error("Unknown label: $label")
    end
end

# Helper: compute convergence rate (slope of log-log plot)
function compute_convergence_rate(xvals, errors)
    logx = log10.(xvals)
    loge = log10.(errors)
    A = hcat(ones(length(logx)), logx)
    coeffs = A \ loge
    return coeffs[2]  # slope
end
# Example usage with cubic_spline
#interpolation_convergence_test(cubic_spline)
#interpolation_convergence_test(newton_interpolation)
interpolation_convergence_test(hermite_interpolation; use_deriv=true)
#interpolation_convergence_test(weno5_interpolate_at)
