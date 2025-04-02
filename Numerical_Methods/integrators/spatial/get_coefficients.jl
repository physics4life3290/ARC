




function get_weights(weights, rule_type, stencil_size)

    # Select integration coefficients
    coefficients = nothing
    
    if weights in [:default, :simpsons]

        if rule_type in [:default, :classic] 

            if stencil_size > 5
                error("The stencil_size of integration must be <5 for Classic coefficients.")
            end

            return simpsons_weights[:classic, stencil_size]

        elseif rule_type == :lownoise 

            if stencil_size < 5
                error("The stencil_size of integration must be ≥5 for Low Noise coefficients.")
            end

            return simpsons_weights[:lownoise, stencil_size]

        end

    elseif weights == :booles

        if rule_type in [:default, :classic] && stencil_size > 7

            if stencil_size > 7
                error("The stencil_size of integration must be <7 for Classic coefficients.")
            end

            return booles_weights[:classic, stencil_size]

        elseif rule_type == :lownoise 
            
            if stencil_size < 7
                error("The stencil_size of integration must be ≥7 for Low Noise coefficients.")
            end

            return booles_weights[:lownoise, stencil_size]

        end

    end

end