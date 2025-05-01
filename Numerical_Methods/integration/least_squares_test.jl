using Plots

# Function to compute least squares polynomial fit
function polynomial_least_squares(x, y, degree)
    X = hcat([x.^i for i in 0:degree]...)  # Proper matrix construction
    β = (X' * X) \ (X' * y)
    return β
end

# Function to evaluate polynomial at a point
function evaluate_polynomial(β, x)
    sum(β[i+1] * x^i for i in 0:length(β)-1)
end

# Example data
x = collect(randn(100))                                
#y = 2 .* (x.+ 0.5 * randn(length(x))) .+ 3 .* (x.+ 0.5 * randn(length(x))).^2 .+ 4 .* (x.+ 0.5 * randn(length(x))).^3   # Noisy quadratic data
y = 5.7 .* (1 .- (x ./ maximum(x)).^2)
degree = 2  # Fitting a degree-2 (quadratic) polynomial
coefficients = polynomial_least_squares(x, y, degree)
println(coefficients)
# Now create fine x points to plot smooth curve
x_fit = range(minimum(x), maximum(x), length=200)
y_fit = [evaluate_polynomial(coefficients, xi) for xi in x_fit]

# Plotting
scatter(x, y, label="Data Points", legend=:topleft)
plot!(x_fit, y_fit, label="Fitted Polynomial", linewidth=2)
