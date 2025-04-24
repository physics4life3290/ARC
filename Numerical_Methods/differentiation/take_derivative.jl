




include("derivatives_include.jl")

function differentiate(x::Vector{Float64}, y::Vector{Float64}; power=1, order=2, difference_type=:default, interp=true)
    
    coefficients = get_derivative_coefficients(difference_type, power, order)

    #println("The coefficients used to compute the integral are: ", coefficients)

    h = (x[end] - x[1]) / (length(x) - 1)
    
    h_optimal = sqrt(eps(Float64)) * abs(sum(x)/length(x))  
    if h < h_optimal  
        @warn "Step size too small, potential numerical issues."
        h = h_optimal
    end

    # Ensure the number of points and coefficients match
    if length(coefficients) != length(x)
        if interp === true
            interp_func(var) = interpolate_to_point(x, y, var)
        end
        interp_points = collect(range(x[1], x[end], length(coefficients)))
        h = (interp_points[end] - interp_points[1]) / (length(coefficients) - 1)
        y = interp_func.(interp_points)  # Ensure broadcasting is applied properly    
    end

    
    # Compute the derivative
    if difference_type in [:default, :forward, :backward]
        return sum(coefficients .* y) /  h
    elseif difference_type == :central
        return sum(coefficients .* y) / (h)
    end
    
end

#=
# Example usage
x = [0.0, 0.01, 0.02, 0.03]  # Forward/backward difference grid
x = x .+ 2.0
f(var) = var^7
println(7 * 2 ^ 6)
y = f.(x)

println(differentiate(x, y, difference_type=:forward, order=4))

println(-1 * differentiate(x, y, difference_type=:backward, order=4))

println(differentiate(x, y, difference_type=:central, order=2))
=#