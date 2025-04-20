

include("interpolation_include.jl")


function interpolate_to_point(ind_var, dep_var, eval_ind_var; 
                              method::Symbol = :div_diff, 
                              mode::Symbol = :minimal, 
                              d_dep_var = nothing)

    # --- Sanity checks ---
    if eval_ind_var < minimum(ind_var) || eval_ind_var > maximum(ind_var)
        @warn "Extrapolating outside the given data range. Results may be unreliable."
    end

    if !issorted(ind_var)
        sortperm_ = sortperm(ind_var)
        ind_var = ind_var[sortperm_]
        dep_var = dep_var[sortperm_]
        if d_dep_var !== nothing
            d_dep_var = d_dep_var[sortperm_]
        end
    end

    if length(unique(ind_var)) < length(ind_var)
        @warn "Duplicate values in independent variable; interpolation may be unstable."
    end

    if length(ind_var) < 3
        error("At least three points are required for interpolation.")
    end

    if length(ind_var) != length(dep_var)
        error("Independent and dependent variable arrays must be the same length.")
    end

    if method == :hermite
        if d_dep_var === nothing
            error("Hermite interpolation requires derivative values.")
        elseif length(ind_var) < 4
            @warn "Hermite interpolation is more stable with at least 4 points."
        end
    end

    # --- Interpolation selection ---
    if mode == :minimal
        return _select_interpolation_method(method, ind_var, dep_var, eval_ind_var, d_dep_var)
    else
        error("Only :minimal mode is implemented.")
    end
end

function _select_interpolation_method(method, x, y, x_eval, dy=nothing)
    if method == :div_diff
        return newton_interpolation(x, y, x_eval)
    elseif method == :spline
        return cubic_spline_interpolation(x, y, x_eval)
    elseif method == :hermite
        return hermite_interpolation(x, y, dy, x_eval)
    elseif method == :WENO
        return weno5_interpolate_at(y, x, x_eval)
    else
        error("Unknown method: $method. Supported: :div_diff, :spline, :hermite, :WENO")
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