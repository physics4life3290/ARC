



function cubic_spline(x::AbstractVector{T}, y::AbstractVector{T}, interp_x; throw_on_bounds::Bool = false) where T <: Real
    n = length(x)
    h = diff(x)

    # Step 1: Compute α for the tridiagonal system
    α = [zero(T); 3.0 * (y[3:end] .- y[2:end-1]) ./ h[2:end] .- 3.0 * (y[2:end-1] .- y[1:end-2]) ./ h[1:end-1]]

    # Step 2: Solve tridiagonal system for c
    l = ones(T, n)
    μ = zeros(T, n)
    z = zeros(T, n)

    for i in 2:n-1
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    c = zeros(T, n)
    b = zeros(T, n - 1)
    d = zeros(T, n - 1)

    for j in n-1:-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end

    # Step 3: Return interpolant function
    function spline_interp(xi::T)
        if throw_on_bounds && (xi < x[1] || xi > x[end])
            error("Interpolation point $xi out of bounds.")
        end

        # Use binary search to find the interval
        i = searchsortedlast(x, xi)
        i = clamp(i, 1, n - 1)

        dx = xi - x[i]
        return y[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3
    end

    return spline_interp(interp_x)
end

#=
x = collect(0.0:0.1:10.0)
y = sin.(x)
spline(var) = cubic_spline(x, y, var)

# Single interpolation
println(spline(1.75))
=#