




function get_derivative_coefficients(difference_type, power, order)
    if difference_type in [:default, :forward]
        coefficients = get(forward_difference_weights, (power, order), nothing)
        if coefficients === nothing
            error("Forward difference not available for power=$power, order=$order")
        end
    elseif difference_type == :backward
        coefficients = get(forward_difference_weights, (power, order), nothing)
        if coefficients === nothing
            error("Backward difference not available for power=$power, order=$order")
        end
        if power == 1
            coefficients = -1 .* coefficients  # Negate for backward difference
        end
    elseif difference_type == :central
        coefficients = get(central_difference_weights, (power, order), nothing)
        if coefficients === nothing
            error("Central difference not available for power=$power, order=$order")
        end
    else
        error("Unknown difference type: $difference_type")
    end
    return coefficients
end