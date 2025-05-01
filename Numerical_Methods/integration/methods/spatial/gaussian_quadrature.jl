




# Same: Legendre polynomial and derivative
function legendre(n, x)
    if n == 0 return (1.0, 0.0) end
    if n == 1 return (x, 1.0) end

    Pnm2, Pnm1 = 1.0, x
    dPnm2, dPnm1 = 0.0, 1.0

    for k = 2:n
        Pn = ((2k - 1) * x * Pnm1 - (k - 1) * Pnm2) / k
        dPn = ((2k - 1) * (Pnm1 + x * dPnm1) - (k - 1) * dPnm2) / k
        Pnm2, Pnm1 = Pnm1, Pn
        dPnm2, dPnm1 = dPnm1, dPn
    end

    return (Pnm1, dPnm1)
end

# Same: find roots
function legendre_roots(n; tol=1e-14, maxiter=100)
    roots = zeros(n)
    for i = 1:div(n + 1, 2)
        x = cos(pi * (i - 0.25) / (n + 0.5))
        for _ = 1:maxiter
            Pn, dPn = legendre(n, x)
            dx = -Pn / dPn
            x += dx
            if abs(dx) < tol break end
        end
        roots[i] = x
        roots[n + 1 - i] = -x
    end
    return roots
end

# Same: compute weights
function legendre_weights(n, roots)
    weights = zeros(n)
    for i in 1:n
        _, dPn = legendre(n, roots[i])
        weights[i] = 2 / ((1 - roots[i]^2) * dPn^2)
    end
    return weights
end

#  Linear interpolation utility
function interpolate(x_vals, y_vals, x)
    for i in 1:length(x_vals)-1
        x0, x1 = x_vals[i], x_vals[i+1]
        if x0 ≤ x ≤ x1
            y0, y1 = y_vals[i], y_vals[i+1]
            return y0 + (x - x0) / (x1 - x0) * (y1 - y0)
        end
    end
    return 0.0  # return 0 if x is outside bounds (simple fallback)
end

#  Gaussian quadrature using sampled function data
function gauss_quad(x_vals, f_vals, n, interp_method=:default)
    a, b = x_vals[1], x_vals[end]
    roots = legendre_roots(n)
    weights = legendre_weights(n, roots)

    total = 0.0
    for i in 1:n
        x_mapped = (b - a) / 2 * roots[i] + (a + b) / 2
        f_mapped = interpolation_dispatch(x_vals, f_vals, x_mapped)
        total += weights[i] * f_mapped
    end
    return (b - a) / 2 * total
end


#=
# -------------------------------------
#  Example with array input
# -------------------------------------
x = range(0, 1, length=100)
f = x .^ 2
approx = gauss_quad(x, f, 5)
println("∫₀¹ x² dx ≈ $approx (exact: 1/3)")
=#