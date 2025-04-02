
#using LinearAlgebra

function cubic_spline_interp(x_vals, y_vals, x)
    n = length(x_vals)
    if n < 3
        throw(ArgumentError("At least three points are required for cubic spline interpolation."))
    end

    # Step 1: Compute step sizes (h)
    h = diff(x_vals)

    # Step 2: Compute the alpha coefficients
    alpha = zeros(n)
    for i in 2:n-1
        alpha[i] = 3 * ((y_vals[i+1] - y_vals[i]) / h[i] - (y_vals[i] - y_vals[i-1]) / h[i-1])
    end

    # Step 3: Solve the tridiagonal system for c-coefficients
    l = ones(n)
    mu = zeros(n)
    z = zeros(n)

    for i in 2:n-1
        l[i] = 2 * (x_vals[i+1] - x_vals[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
    end

    # Enforce natural spline conditions (c[1] = c[n] = 0)
    c = zeros(n)
    b = zeros(n - 1)
    d = zeros(n - 1)

    for j in (n-1):-1:1
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (y_vals[j+1] - y_vals[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end

    # Step 4: Find the interval for x and interpolate
    idx = searchsortedlast(x_vals, x)
    idx = clamp(idx, 1, n-1)  # Ensure valid indexing

    dx = x - x_vals[idx]
    return y_vals[idx] + b[idx] * dx + c[idx] * dx^2 + d[idx] * dx^3
end

# Minmod limiter function
function minmod(a, b)
    return sign(a) * min(abs(a), abs(b))
end

# Apply limiter to second derivative (c) or first derivative (b)
function limit_spline_derivatives(x_vals, y_vals, b, c, h)
    n = length(x_vals)

    # Apply minmod limiter to the slopes (first derivatives)
    for i in 2:n-1
        slope_left = (y_vals[i] - y_vals[i-1]) / h[i-1]
        slope_right = (y_vals[i+1] - y_vals[i]) / h[i]
        
        # Apply limiter to slope
        limited_slope = minmod(slope_left, slope_right)

        # Correct the second derivative based on the limited slope
        if limited_slope != 0
            c[i] = (limited_slope - (y_vals[i+1] - y_vals[i-1]) / (h[i-1] + h[i])) / 2
        end
    end

    return b, c
end

function cubic_spline_interp_with_limiters(x_vals, y_vals, x)
    n = length(x_vals)
    if n < 3
        throw(ArgumentError("At least three points are required for cubic spline interpolation."))
    end

    # Step 1: Calculate step sizes (h)
    h = diff(x_vals)

    # Step 2: Compute coefficients for the tridiagonal system
    alpha = zeros(n - 1)
    for i in 2:n-1
        alpha[i] = 3 * ((y_vals[i+1] - y_vals[i]) / h[i] - (y_vals[i] - y_vals[i-1]) / h[i-1])
    end

    # Step 3: Solve the tridiagonal system for second derivatives (c)
    l = ones(n)
    mu = zeros(n)
    z = zeros(n)

    for i in 2:n-1
        l[i] = 2 * (x_vals[i+1] - x_vals[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
    end

    c = zeros(n)
    b = zeros(n - 1)
    d = zeros(n - 1)

    for j in (n-1):-1:1
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (y_vals[j+1] - y_vals[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end

    # Step 4: Apply the limiter to prevent steepening and flattening
    b, c = limit_spline_derivatives(x_vals, y_vals, b, c, h)

    # Step 5: Find the interval for x and interpolate
    idx = searchsortedlast(x_vals, x)
    if idx == length(x_vals)  # Handle boundary cases
        return y_vals[end]
    elseif idx == 0
        return y_vals[1]
    else
        dx = x - x_vals[idx]
        return y_vals[idx] + b[idx] * dx + c[idx] * dx^2 + d[idx] * dx^3
    end
end


using LinearAlgebra

# Function to compute the B-spline basis function N(i, k, t) recursively
function bspline_basis(i, k, t, knots)
    if k == 0
        return knots[i] ≤ t < knots[i+1] ? 1.0 : 0.0
    else
        denom1 = knots[i+k] - knots[i]
        denom2 = knots[i+k+1] - knots[i+1]
        
        term1 = denom1 > 0 ? ((t - knots[i]) / denom1) * bspline_basis(i, k-1, t, knots) : 0.0
        term2 = denom2 > 0 ? ((knots[i+k+1] - t) / denom2) * bspline_basis(i+1, k-1, t, knots) : 0.0

        return term1 + term2
    end
end

# Function to evaluate a cubic B-spline curve at a given t
function cubic_bspline(t, control_points, knots)
    n = length(control_points)
    k = 3  # Cubic B-spline (degree 3)
    
    result = zero(control_points[1])  # Initialize to correct dimensionality
    for i in 1:n
        coeff = bspline_basis(i, k, t, knots)
        result += coeff * control_points[i]
    end
    
    return result
end

#=
# Example Usage:
control_points = [0.0, 1.0, 2.0, 3.0, 4.0]  # Example 1D control points
knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]  # Uniform knot vector

t_values = range(0.0, stop=3.0, length=100)
b_spline_values = [cubic_bspline(t, control_points, knots) for t in t_values]

using Plots
println("B-Spline evaluated at sampled points:", b_spline_values)
plot(t_values, b_spline_values)
=#