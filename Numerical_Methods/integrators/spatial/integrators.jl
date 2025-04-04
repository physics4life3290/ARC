





include("integrate_include.jl")


function integrate(ind_var, integrand; weights=:default, rule_type=:default, stencil_size=3, interp_method::Symbol=:default)
    println("The points coming into the integral are: ", ind_var)
    println("The values of the integrand at these points are: ", integrand)
    println("The weights being used to compute this integral are from the ", weights, " method.")
    println("The class of coefficients being used are: ", rule_type)
    println("The stencil size is: ", stencil_size)
    
    if rule_type != :quadrature
        coefficients =  get_weights(weights, rule_type, stencil_size)
        println("The coefficients being used to compute this integral are: ", coefficients)
        # Ensure coefficients exist
        if coefficients === nothing
            error("Unsupported rule type or stencil_size: ($rule_type, $stencil_size)")
        end
        return calculate_integral(ind_var, integrand, coefficients, interp_method=interp_method)
    end
    
end


#=
ω = 1.0
f(var) = sin(ω * var)
x = range(0, π/2, length=3) |> collect
y = f.(x)

result = integrate(x, y; weights=:simpsons, rule_type=:classic, stencil_size=4)
#println(result)
#println("Error from answer is: ", (result[2]-0.45970)/0.45970)
=#