using Plots
include("../interpolation_include.jl")

# Define test functions and their derivatives
f_smooth(var) = var^3
df_smooth(var) = 3 * var^2
f_oscillatory(var) = sin(var)
df_oscillatory(var) = cos(var)
f_square_wave(var) = sign(sin(var))
df_square_wave(var) = cos(var) * (var > 0 ? 1 : -1)

# Define independent variable
independent_var = collect(range(0.0, stop=10.0, length=100))

# Compute dependent variables
smooth_dep_var = f_smooth.(independent_var)
dsmooth_dep_var = df_smooth.(independent_var)
oscillatory_dep_var = f_oscillatory.(independent_var)
doscillatory_dep_var = df_oscillatory.(independent_var)
square_wave_dep_var = f_square_wave.(independent_var)
dsquare_wave_dep_var = df_square_wave.(independent_var)

# Store variable data and labels
test_variables = Dict(
    "smooth" => (smooth_dep_var, dsmooth_dep_var),
    "oscillatory" => (oscillatory_dep_var, doscillatory_dep_var),
    "square_wave" => (square_wave_dep_var, dsquare_wave_dep_var),
)


lengths = [10, 20, 50, 100]

function interpolation_convergence_test(interp_func::Function; use_deriv::Bool=false)
    all_results = Dict{String, Dict{Int, Vector{Float64}}}()
    interp_xs = Dict{Int, Vector{Float64}}()

    for (label, (dep_var, ddep_var)) in test_variables
        all_results[label] = Dict{Int, Vector{Float64}}()

        for len in lengths
            interp_ind_var = collect(range(0.5, stop=9.5, length=len))
            interp_dep_var = Vector{Float64}(undef, len)

            for i in 1:len
                x_eval = interp_ind_var[i]

                interp_dep_var[i] = use_deriv ?
                    interp_func(independent_var, dep_var, ddep_var, x_eval) :
                    interp_func(independent_var, dep_var, x_eval)
            end

            all_results[label][len] = interp_dep_var
            interp_xs[len] = interp_ind_var
        end
    end

    # Plotting
    for (label, result_dict) in all_results
        plt = plot(title="Interpolations for $label function", xlabel="x", ylabel="Interpolated f(x)")
        for len in lengths
            plot!(plt, interp_xs[len], result_dict[len], label="n=$len", lw=2)
        end
        display(plt)
        # savefig(plt, "interpolation_$(label).png")
    end
end

# Example usage with cubic_spline
#interpolation_convergence_test(cubic_spline)
#interpolation_convergence_test(newton_interpolation)
#interpolation_convergence_test(hermite_interpolation; use_deriv=true)
interpolation_convergence_test(weno5_interpolate_at)
