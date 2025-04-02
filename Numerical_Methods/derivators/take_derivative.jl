




include("derivatives_include.jl")

function differentiate(x, y; power=1, order=2, difference_type=:default)
    
    h = (x[end] - x[1]) / (length(x) - 1)
    
    coefficients = get_derivative_coefficients(difference_type, power, order)

    # Ensure the number of points and coefficients match
    if length(coefficients) != length(x)
        interp_points = collect(range(x[1], x[end], length(coefficients)))
        h = (interp_points[end] - interp_points[1]) / (length(coefficients) - 1)
        y_interp(var) = cubic_spline_interp(x, y, var)
        y = y_interp.(interp_points)
    end

    # Compute the derivative
    if difference_type in [:default, :forward, :backward]
        return dot(coefficients, y) / h
    elseif difference_type == :central
        return dot(coefficients, y) / (2 * h)
    end
end



#=
# Example usage
x = [0.0, 0.33333,0.66667, 1]  # Forward/backward difference grid
f(var) = var^3
y = f.(x)

println(differentiate(x, y, difference_type=:forward, order=6))

x = [0.0, 0.5, 1.0, 1.5]  # Four points needed for backward difference (order=3)
f(var) = var^3
y = f.(x)

println(differentiate(x, y, difference_type=:backward, order=3))
=#