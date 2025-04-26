function weno5_interpolate_at(x::Vector{Float64}, f::Vector{Float64}, x_interp::Float64; interpolation_type::Symbol=:left)
    N = length(f)
    dx = x[2] - x[1]  # Assumes uniform spacing

    # Find cell index i such that x[i] ≤ x_interp < x[i+1]
    idx = findfirst(z -> z > x_interp, x)
    if idx === nothing
        error("x_interp = $x_interp is out of the domain of x.")
    end
    i = idx - 1
    if i < 3 || i > N - 2
        error("x_interp = $x_interp is too close to the boundaries for WENO5 interpolation.")
    end

    # Relative position in the cell [x[i], x[i+1]] (unused in this simple WENO5)
    # ξ = (x_interp - x[i]) / dx

    # Load stencil values
    f0, f1, f2, f3, f4 = f[i-2], f[i-1], f[i], f[i+1], f[i+2]

    ϵ = 1e-6  # small number to avoid division by zero

    # Candidate interpolants (left-biased)
    if interpolation_type == :left
        q0 = (2*f0 - 7*f1 + 11*f2) / 6
        q1 = (-f1 + 5*f2 + 2*f3) / 6
        q2 = (2*f2 + 5*f3 - f4) / 6
    elseif interpolation_type == :centered
        # A simple centered approximation across the stencil
        q0 = (f0 + f1 + f2) / 3
        q1 = (f1 + f2 + f3) / 3
        q2 = (f2 + f3 + f4) / 3
    else
        error("Invalid interpolation type. Use :left or :centered.")
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

    # Interpolated value
    return w0*q0 + w1*q1 + w2*q2
end

#=
# Example usage
x = collect(LinRange(-2.0, 2.0, 1000))
f = [xi < 0 ? 1.0 : xi^2 for xi in x]

for interp_type in (:left, :centered)
    y = weno5_interpolate_at(x, f, 0.5; interpolation_type=interp_type)
    println("Interpolation (", interp_type, ") at x=0.5: ", y)
end
=#