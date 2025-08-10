




# Minmod slope limiter
function minmod(a, b)
    if a * b <= 0
        return 0.0
    else
        return sign(a) * min(abs(a), abs(b))
    end
end

function vanleer(a, b)
    if a * b <= 0
        return 0.0
    else
        return (2 * a * b) / (a + b)
    end
end

function superbee(a, b)
    if a * b <= 0
        return 0.0
    else
        s = sign(a)
        return s * max(0.0, min(2 * abs(a), abs(b)), min(abs(a), 2 * abs(b)))
    end
end


# Compute slopes for primitive variables with minmod limiter
function compute_slopes(var)
    n = length(var)
    slopes = zeros(n)
    for i in 2:n-1
        dl = var[i] - var[i-1]
        dr = var[i+1] - var[i]
#        slopes[i] = minmod(dl, dr)
#        slopes[i] = superbee(dl, dr)
        slopes[i] = vanleer(dl, dr)
    end
    slopes[1] = 0.0
    slopes[end] = 0.0
    return slopes
end

# Reconstruction at interfaces: left and right states
function reconstruct(var, slopes)

    n = length(var)
    # Left state at interface i+1/2 comes from cell i
    var_L = zeros(n+1)
    # Right state at interface i+1/2 comes from cell i+1
    var_R = zeros(n+1)

    # Interior interfaces
    for i in 2:n
        var_L[i] = var[i] + 0.5 * slopes[i]
        var_R[i] = var[i] - 0.5 * slopes[i]
    end

    # Boundary interfaces (use cell values, no extrapolation)
    var_L[1] = var[1]
    var_R[1] = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end

function reconstruct_parabolic(var)
    n = length(var)

    # Left and Right reconstructed states at interfaces (n+1 interfaces)
    var_L = zeros(n+1)
    var_R = zeros(n+1)

    # Compute parabolic interface values using a 3-point stencil
    for i in 2:n-1
        # Fit a parabola through var[i-1], var[i], var[i+1]
        dq_minus = var[i] - var[i-1]
        dq_plus = var[i+1] - var[i]

        # Centered difference with limiter (basic minmod to reduce overshoot)
        slope = minmod(dq_minus, dq_plus)

        # Compute curvature
        d2q = var[i+1] - 2*var[i] + var[i-1]

        # Parabolic reconstruction
        var_L[i+1] = var[i] + 0.5 * slope - (1/6) * d2q
        var_R[i]   = var[i] - 0.5 * slope - (1/6) * d2q
    end

    # Boundary extrapolation (copy cell-centered value)
    var_L[1]   = var[1]
    var_R[1]   = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end


function steepening(ρ, p, slope, iters)
    for i in 2:iters-1
        dCρ = 0.5 * (ρ[i+1] - ρ[i-1])
        ΔpL = abs(p[i] - p[i-1])
        ΔpR = abs(p[i+1] - p[i])
        ΔρL = abs(ρ[i]-ρ[i-1])
        ΔρR = abs(ρ[i+1]-ρ[i])
        steepening_coef = 0.0

        if (ΔρL > 0.02 * ρ[i]) && (ΔpL < 0.1 * p[i]) &&
            (ΔρR > 0.02 * ρ[i]) && (ΔpR < 0.1 * p[i])
            steepening_coef = 0.3  # strengthen slope slightly
        end
        slope[i] += steepening_coef * dCρ
    end
    
    return slope
end

function MUSCL_Step!(UserInput, _grid, W, U, dt)
    
    nx = _grid.xcoord.total_zones
    γ = UserInput.secondary_input.γ
    dx = _grid.xcoord.spacing

    # Compute slopes for primitive variables
    slope_ρ = compute_slopes(W.density_centers)
    slope_u = compute_slopes(W.velocity_centers)
    slope_p = compute_slopes(W.pressure_centers)
    for i in 2:nx-1
        shock_strength = abs(W.density_centers[i+1] - 2*W.density_centers[i] + W.density_centers[i-1])  # curvature of density
        flattening_coef = clamp(1.0 - 10.0 * shock_strength, 0.0, 1.0)
        slope_ρ[i] = slope_ρ[i] * flattening_coef
    end

    slope_ρ = steepening(W.density_centers, W.pressure_centers, slope_ρ, _grid.xcoord.total_zones)
    
    # Reconstruct left and right states at interfaces
    ρL, ρR = reconstruct(W.density_centers, slope_ρ)
    uL, uR = reconstruct(W.velocity_centers, slope_u)
    pL, pR = reconstruct(W.pressure_centers, slope_p)
    #ρL, ρR = reconstruct_parabolic(W.density_centers)
    #uL, uR = reconstruct_parabolic(W.velocity_centers)
    #pL, pR = reconstruct_parabolic(W.pressure_centers)

    apply_boundary_conditions(UserInput, U, _grid)
    # Fluxes at interfaces
    F1 = zeros(nx+1)
    F2 = zeros(nx+1)
    F3 = zeros(nx+1)
    Threads.@threads for i in 2:nx
        @inbounds begin
            f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
            #f1, f2, f3 = ExactRiemannSolve!((ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i]),γ)
            F1[i] = f1
            F2[i] = f2
            F3[i] = f3

        end
    end

    # Update conserved variables (interior cells)
    Threads.@threads for i in 2:nx-1
        @inbounds begin
            U.density_centers[i] -= dt/dx * (F1[i+1] - F1[i])
            U.momentum_centers[i] -= dt/dx * (F2[i+1] - F2[i])
            U.total_energy_centers[i] -= dt/dx * (F3[i+1] - F3[i])

        end
    end
end