




include("derivatives_include.jl")

function differentiate(x::Vector{Float64}, y::Vector{Float64}; power=1, order=2, difference_type=:default, interp=true)
    
    #println("The points coming into the function are: ", x)
    #println("The values associated with those points are: ", y)
    #println("You are taking a derivative or order: ", power)
    #println("The order of accuracy of this derivative is: ", order)
    #println("The option to interpolate is set to: ", interp)

    coefficients = get_derivative_coefficients(difference_type, power, order)

    #println("The coefficients used to compute the integral are: ", coefficients)

    h = (x[end] - x[1]) / (length(x) - 1)
    
    h_optimal = sqrt(eps(Float64)) * abs(sum(x)/length(x))  
    if h < h_optimal  
        @warn "Step size too small, potential numerical issues."
        h = h_optimal
    end

    #println("The initial step is: ", h)
    #println("The machines optimal step is: ", h_optimal)
    
    #println("The conditional for interpolation is: ", length(coefficients) != length(x))
    # Ensure the number of points and coefficients match
    if length(coefficients) != length(x)
        if interp === true
            interp_func(var) = interpolate_to_point(x, y, var)
        end
        interp_points = collect(range(x[1], x[end], length(coefficients)))
        #println("The points used to calculate the derivative are: ", interp_points)
        h = (interp_points[end] - interp_points[1]) / (length(coefficients) - 1)
        #println("The new step size for the computed derivative is: ", h)
        y = interp_func.(interp_points)  # Ensure broadcasting is applied properly
        #println("The new y values used for differentiation: ", y)    
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