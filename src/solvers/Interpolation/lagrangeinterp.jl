









function lagrange_interp(xdata::Vector{Float64}, ydata::Vector{Float64}, x::Float64;
                         degree::Int=2, side::Symbol=:center)
    N = length(xdata)
    npts = degree + 1
    if npts > N
        error("Degree too high for number of data points.")
    end

    # Find the index nearest to x
    idx = findfirst(xi -> xi ≥ x, xdata)
    if idx === nothing
        idx = N
    end

    # Choose indices of interpolation points
    if side == :left
        istart = max(1, idx - npts + 1)
        iend = istart + npts - 1
    elseif side == :right
        iend = min(N, idx + npts - 1)
        istart = iend - npts + 1
    else # center
        half = npts ÷ 2
        istart = max(1, idx - half)
        iend = istart + npts - 1
        if iend > N
            iend = N
            istart = N - npts + 1
        end
    end

    xs = xdata[istart:iend]
    ys = ydata[istart:iend]

    # Compute Lagrange polynomial
    Lx = 0.0
    iters = length(xs)
    for i in 1:iters
        li = 1.0
        for j in 1:iters
            if i != j
                li *= (x - xs[j]) / (xs[i] - xs[j])
            end
        end
        Lx += li * ys[i]
    end

    return Lx
end

#=
# Example usage
xdata = 0.0:1.0:5.0
ydata = xdata.^3
x = 2.5

println(lagrange_interp(collect(xdata), ydata, x, degree=3, side=:center))
=#