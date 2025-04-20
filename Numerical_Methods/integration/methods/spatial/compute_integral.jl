




function calculate_integral(bounds, integrand, coefficients; interp_method=:default)
    
    if length(bounds) == length(coefficients)

        h = (bounds[end] - bounds[1]) / (length(bounds) - 1)
        integral_array = h .* coefficients .* integrand
        integral_sum = sum(integral_array) 
    
    elseif length(bounds) != length(coefficients)
        
        if interp_method == :default
            interp_method = :div_diff
        end

        interp_func(var) = interpolate_to_point(bounds, integrand, var, method=interp_method)
        interp_points = collect(range(bounds[1], bounds[end], length(coefficients)))
        interp_y = interp_func.(interp_points)
        # Compute integral
        h = (bounds[end] - bounds[1]) / (length(interp_points) - 1)
        integral_array = h .* coefficients .* interp_y
        integral_sum = sum(integral_array)

    end

    return integral_sum
    
end