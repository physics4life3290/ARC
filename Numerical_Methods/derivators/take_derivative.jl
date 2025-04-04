




include("derivatives_include.jl")

function differentiate(x::Vector{Float64}, y::Vector{Float64}; power=1, order=2, difference_type=:default, interp_func=nothing)

    coefficients = get_derivative_coefficients(difference_type, power, order)

    h = (x[end] - x[1]) / (length(x) - 1)
    
    h_optimal = sqrt(eps(Float64)) * abs(sum(x)/length(x))  
    if h < h_optimal  
        @warn "Step size too small, potential numerical issues."
        h = h_optimal
    end


    # Ensure the number of points and coefficients match
    if length(coefficients) != length(x)
        if interp_func === nothing
            interp_func(ind, dep, var) =  
        interp_points = collect(range(x[1], x[end], length(coefficients)))
        h = (interp_points[end] - interp_points[1]) / (length(coefficients) - 1)
        y_interp(var) = interp_func(x, y, var)
        y = y_interp.(interp_points)
    end

    # Compute the derivative
    if difference_type in [:default, :forward, :backward]
        return dot(coefficients, y) /  order * h
    elseif difference_type == :central
        return dot(coefficients, y) / (2 * h)
    end
end




# Example usage
x = [0.0, π/2, π]  # Forward/backward difference grid
f(var) = sin(var)
y = f.(x)

println(differentiate(x, y, difference_type=:forward, order=4))

println(differentiate(x, y, difference_type=:backward, order=4))


println(differentiate(x, y, difference_type=:central, order=4))