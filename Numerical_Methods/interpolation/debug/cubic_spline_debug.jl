
function spline_interp(xi, x, y, b, c, d)
    i = searchsortedlast(x, xi)
    i = clamp(i, 1, length(x) - 1)

    dx = xi - x[i]
    return y[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3
end

function cubic_spline_debug(x::AbstractVector{T}, y::AbstractVector{T}, interp_ind_var;
                            throw_on_bounds::Bool = false,
                            logfile::String = "Numerical_Methods/interpolation/debug/output/cubic_spline_debug.txt") where T <: Real

    n = length(x)
    h = diff(x)
    α = [zero(T); 3.0 * (y[3:end] .- y[2:end-1]) ./ h[2:end] .- 3.0 * (y[2:end-1] .- y[1:end-2]) ./ h[1:end-1]]

    l = ones(T, n)
    μ = zeros(T, n)
    z = zeros(T, n)

    c = zeros(T, n)
    b = zeros(T, n - 1)
    d = zeros(T, n - 1)

    open(logfile, "w") do io
        try
            write(io, "Cubic Spline Debug Log\n")
            write(io, "======================\n")
            write(io, "n = $n\n")
            write(io, "x = $(x)\n")
            write(io, "y = $(y)\n")
            write(io, "h = $(h)\n")
            write(io, "α = $(α)\n")
            write(io, "Initial l = $(l)\n")
            write(io, "Initial μ = $(μ)\n")
            write(io, "Initial z = $(z)\n")

            for i in 2:n-1
                l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
                μ[i] = h[i] / l[i]
                z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
            end

            write(io, "Final l = $(l)\n")
            write(io, "Final μ = $(μ)\n")
            write(io, "Final z = $(z)\n")

            for j in n-1:-1:1
                c[j] = z[j] - μ[j] * c[j+1]
                b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
                d[j] = (c[j+1] - c[j]) / (3 * h[j])
            end

            write(io, "c = $(c)\n")
            write(io, "b = $(b)\n")
            write(io, "d = $(d)\n")
            write(io, "throw_on_bounds = $throw_on_bounds\n")

        catch e
            write(io, "\n--- ERROR ENCOUNTERED ---\n")
            write(io, "Error Type: $(typeof(e))\n")
            write(io, "Error Message: $(e.msg)\n")
            write(io, "Stacktrace:\n$(catch_backtrace())\n")
            rethrow()
        end
    end

    return spline_interp(interp_ind_var, x, y, b, c, d)
end
