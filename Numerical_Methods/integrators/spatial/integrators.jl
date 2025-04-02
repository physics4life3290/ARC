





include("integrate_include.jl")


function integrate(bounds, integrand; weights=:default, rule_type=:default, stencil_size=3, interp_function=nothing)
    if rule_type != :quadrature
        coefficients =  get_weights(weights, rule_type, stencil_size)
        
        # Ensure coefficients exist
        if coefficients === nothing
            error("Unsupported rule type or stencil_size: ($rule_type, $stencil_size)")
        end

        if interp_function === nothing
            return calculate_integral(bounds, integrand, coefficients)
        elseif interp_function !== nothing
            return calculate_integral(bounds, integrand, coefficients, interp_func=interp_function)
        end
    elseif rule_type == :quadrature
        # change of variable to have bounds be from [-1.1]
        maxbound = maximum(bounds)
        minbound = minimum(bounds)
        bounds = ((2 .* bounds) .- (maxbound + minbound)) ./ (maxbound - minbound)
        println(bounds)

        #   i) need to find legendre polynomials
        #  ii) need to find roots of legendre polynomials, this gives us the nodes that we integrate at
        # iii) need to find derivative of legendre polynomial
        function legendre(n, x)
            if n == 0
                return 1
            elseif n == 1
                return x
            else
                P0 = 1
                P1 = x
                for k in 2:n
                    P2 = ((2*k - 1) * x * P1 - (k - 1) * P0) / k
                    P0, P1 = P1, P2
                end
                return P1
            end
        end
        
        # Example: Calculate P_n(x) for n = 0, 1, 2, and x = 0.5
        x = 0.5
        for n in 0:3
            println("P_$n($x) = ", legendre(n, x))
        end
        
    end
end

#=
ω = 1.0
f(var) = sin(ω * var)
x = range(0, 1, length=3) |> collect
y = f.(x)

result = integrate(x, y; weights=:booles, rule_type=:quadrature, stencil_size=8)
println(result)
println("Error from answer is: ", (result[2]-0.45970)/0.45970)
=#