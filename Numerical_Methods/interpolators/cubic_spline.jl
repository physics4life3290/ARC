





function cubic_spline_interpolation(x::Vector{T}, y::Vector{T}, x_interp::T) where T
    n = length(x) - 1
    h = diff(x)
    alpha = [0.0; 3.0 ./ h[2:end] .* (y[3:end] .- y[2:end-1]) .- 3.0 ./ h[1:end-1] .* (y[2:end-1] .- y[1:end-2])]

    # Solve tridiagonal system for c coefficients
    l = ones(T, n + 1)
    mu = zeros(T, n)
    z = zeros(T, n + 1)
    
    for i in 2:n
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
    end

    c = zeros(T, n + 1)
    b = zeros(T, n)
    d = zeros(T, n)
    
    for j in n:-1:1
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end

    function spline_eval(x_val)
        idx = findlast(x .<= x_val)  # Find the interval
        if idx === nothing || idx == length(x)  # Out of bounds handling
            error("Interpolation point out of range.")
        end
        dx = x_val - x[idx]
        return y[idx] + b[idx] * dx + c[idx] * dx^2 + d[idx] * dx^3
    end

    return spline_eval(x_interp)
end
