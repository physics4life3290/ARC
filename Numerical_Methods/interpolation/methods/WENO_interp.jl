




function weno5_interpolate_at(x::Vector{Float64}, f::Vector{Float64}, x_interp::Float64, interpolation_type::Symbol=:left)
    N = length(f)
    dx = x[2] - x[1]  # Assumes uniform spacing

    # Determine index i such that x[i] ≤ x_interp < x[i+1]
    i = findfirst(z -> z > x_interp, x)
    if i === nothing || i < 3 || i > N - 2
        error("x_interp = $x_interp is out of bounds for WENO5 interpolation.")
    end
    #i -= 1

    # Relative position of x_interp in [x[i], x[i+1]]
    ξ = (x_interp - x[i]) / dx

    # Load stencil
    f0, f1, f2, f3, f4 = f[i-2], f[i-1], f[i], f[i+1], f[i+2]

    ϵ = 1e-6  # to avoid division by zero

    # Candidate interpolants (left-biased or centered)
    if interpolation_type == :left
        q0 = (2*f0 - 7*f1 + 11*f2) / 6
        q1 = (-f1 + 5*f2 + 2*f3) / 6
        q2 = (2*f2 + 5*f3 - f4) / 6
    elseif interpolation_type == :centered
        # Centered version: use symmetric combinations
        q0 = (f0 + f1 + f2) / 3
        q1 = (f1 + f2 + f3) / 3
        q2 = (f2 + f3 + f4) / 3
    else
        error("Invalid interpolation type. Choose either :left or :centered.")
    end

    # Smoothness indicators
    β0 = (13/12)*(f0 - 2*f1 + f2)^2 + (1/4)*(f0 - 4*f1 + 3*f2)^2
    β1 = (13/12)*(f1 - 2*f2 + f3)^2 + (1/4)*(f1 - f3)^2
    β2 = (13/12)*(f2 - 2*f3 + f4)^2 + (1/4)*(3*f2 - 4*f3 + f4)^2

    # Linear weights
    γ0, γ1, γ2 = 0.1, 0.6, 0.3

    # Nonlinear weights
    α0 = γ0 / (ϵ + β0)^2
    α1 = γ1 / (ϵ + β1)^2
    α2 = γ2 / (ϵ + β2)^2
    αsum = α0 + α1 + α2

    w0 = α0 / αsum
    w1 = α1 / αsum
    w2 = α2 / αsum

    # Return weighted combination
    return w0*q0 + w1*q1 + w2*q2
end


#=
# Sample function with a discontinuity
x = collect(LinRange(-π, π, 10000))
f = [xi < 0 ? 1.0 : xi^2 for xi in x]

x_interp = 0.5
y_interp = weno5_interpolate_at(f, x, x_interp, :left)
println("Left-biased Interpolated value at x = $x_interp: $y_interp")

x_interp = 0.5
y_interp = weno5_interpolate_at(f, x, x_interp, :centered)
println("Centered Interpolated value at x = $x_interp: $y_interp")
=#