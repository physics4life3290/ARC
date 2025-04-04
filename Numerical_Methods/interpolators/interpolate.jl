

include("interpolation_include.jl")


function interpolate_to_point(ind_var, dep_var, eval_ind_var; method::Symbol=:div_diff, mode=:minimal, d_dep_var=nothing)

        # Check for extrapolation
        if eval_ind_var < minimum(ind_var) || eval_ind_var > maximum(ind_var)
            @warn "Extrapolating outside the given data range. Results may be unreliable."
        end

        # Check for values in ascending order within the array
        if issorted(ind_var) != true
            ind_var = sort!(ind_var)
        end

        # Check for duplicate points
        if length(ind_var) < length(unique(ind_var))
            @warn "Duplicate points detected in the independent variable... \nResults may be incorrect!"
        end

        if length(ind_var) < 3 
            error("At least three points require when interpolating...")
        end

        if length(ind_var) != length(dep_var)
            error("Independent Variable and Dependent Variable must have same number of points...")
        end

        if method == :hermite && d_dep_var === nothing
            error("You must provide values for the derivative...")
        end

        if method == :hermite && length(ind_var) < 4
            @warn "We recommend at least 4 points for interpolation when using the Hermite method..."
        end

        if mode == :minimal
            if method == :div_diff
                return  newton_interpolation(ind_var, dep_var, eval_ind_var)
            elseif method == :spline
                return cubic_spline_interpolation(ind_var, dep_var, eval_ind_var)
            elseif method == :hermite
                return hermite_interpolation(ind_var, dep_var, d_dep_var, eval_ind_var)
            end
        end
end


#= Example usage:
x = [1.0, 2.0, 3.0, 4.0]
y = [1.0, 8.0, 27.0, 64.0] # y = x^2
dy = [3.0, 12.0, 27.0, 48.0]

x_eval = 1.5
interp_value = newton_interpolation(x, y, x_eval)
error_estimate = newton_interpolation_error(x, y, x_eval)

println("Interpolated value at x=$x_eval: ", interp_value)
println("Estimated error at x=$x_eval: ", error_estimate)
=#
#=
# Example usage:
x = [0.0, 1.0, 2.0, 3.0, 4.0]
y = [0.0, 1.0, 8.0, 27.0, 64.0]
x_interp = 1.5
y_interp = cubic_spline_interpolation(x, y, x_interp)
println(y_interp)



# Example Usage
x_vals = [0.0, 1.0, 2.0]      # Given x-values
f_vals = [0.0, 1.0, 8.0]      # f(x) = x^3
df_vals = [0.0, 3.0, 12.0]    # f'(x) = 3x^2

x_interp = 2.5
result = hermite_interpolation(x_vals, f_vals, df_vals, x_interp)

println("Interpolated value at x=$x_interp: ", result)
println("Expected x^3 value at x=$x_interp: ", x_interp^3)


println(interpolate_to_point(x_vals, f_vals, 2.5, method=:hermite, d_dep_var=df_vals))
=#