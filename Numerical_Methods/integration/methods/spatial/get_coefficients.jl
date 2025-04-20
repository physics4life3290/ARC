




function get_weights(weights, rule_type, stencil_size)
    println("The weights we want to use are: ", weights)
    println("We want them to be of the: ", rule_type, " variety")
    println("The stencil size is: ", stencil_size)
    # Select integration coefficients
    coefficients = nothing
    
    if weights in [:default, :simpsons]

        if rule_type in [:default, :classic] 

            if stencil_size > 5
                error("For integration using Simpson's rule, with classic coefficients, the stencil size should be 3 or 4...")
            end

            return simpsons_weights[:classic, stencil_size]

        elseif rule_type == :lownoise 

            if stencil_size < 5
                error("For integration using Simpson's rule, with low noise coefficients, the stencil size should be 5, 6, or 7...")
            end

            return simpsons_weights[:lownoise, stencil_size]

        end

    elseif weights == :booles
        
        if rule_type in [:default, :classic] 
            
            if stencil_size > 7
                error("For integration using Boole's rule, with classic coefficients, the stencil size should be 5 or 6...")
            end

            return booles_weights[:classic, stencil_size]

        elseif rule_type == :lownoise 
            
            if stencil_size < 7
                error("For integration using Boole's rule, with low noise coefficients, the stencil size should be 7, 8, or, 9...")
            end

            return booles_weights[:lownoise, stencil_size]

        end

    end

end