





function newton_divided_differences(x, y)
    n = length(x)
    F = zeros(n, n)
    F[:, 1] = y  # First column is just the y-values
    
    for j in 2:n
        for i in 1:(n-j+1)
            F[i, j] = (F[i+1, j-1] - F[i, j-1]) / (x[i+j-1] - x[i])
        end
    end
    return F[1, :], F
end

function newton_interpolation(x, y, x_eval)

    coeffs, _ = newton_divided_differences(x, y)
    n = length(coeffs)
    
    # Evaluate the Newton polynomial at x_eval
    p = coeffs[n]
    for k in (n-1):-1:1
        p = coeffs[k] + (x_eval - x[k]) * p
    end
    return p
end

function newton_interpolation_error(x, y, x_eval)
    n = length(x)
    _, F = newton_divided_differences(x, y)
    
    # Estimate the (n+1)th derivative using the last divided difference
    max_derivative = abs(F[1, n])
    
    # Compute the error bound
    error_term = 1.0
    for xi in x
        error_term *= (x_eval - xi)
    end
    
    error_bound = abs(error_term) * max_derivative / factorial(n)
    return error_bound
end

#=
# Define function to interpolate
f(x) = sin(x)

# Sample data points
x_nodes = [0.0, π/6, π/4, π/3]
y_nodes = f.(x_nodes)

# Point to evaluate the interpolation
x_eval = π/5  # ≈ 0.628

# Compute interpolated value
p_val = newton_interpolation(x_nodes, y_nodes, x_eval)

# Compute true value
true_val = f(x_eval)

# Compute error estimate
err_est = newton_interpolation_error(x_nodes, y_nodes, x_eval)

# Output results
println("Interpolating sin(x) at x = $x_eval")
println("Interpolated value     : $p_val")
println("True value             : $true_val")
println("Absolute error         : $(abs(p_val - true_val))")
println("Estimated error bound  : $err_est")
=#