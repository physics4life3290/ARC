




function calculate_integral(bounds, integrand, coefficients; interp_func=nothing)

    if length(bounds) == length(coefficients)
        h = (bounds[end] - bounds[1]) / (length(bounds) - 1)
        integral_array = h .* coefficients .* integrand
        integral_sum = sum(integral_array) 
    elseif length(bounds) != length(coefficients)
        interp_points = collect(range(bounds[1], bounds[end], length(coefficients)))
        interp_y = interp_func.(interp_points)
        # Compute integral
        h = (bounds[end] - bounds[1]) / (length(interp_points) - 1)
        integral_array = h .* coefficients .* interp_y
        integral_sum = sum(integral_array)
    end

    return integral_sum
end