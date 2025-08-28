


# I am updating every ounce of code to support generalization for anyone to use. #

# So, serial, threads, and GPU modes should have boolean features like verbose, debug, and benchmark 

function flattening(W, slope_ρ)
    nx = length(W.density_centers)
    for i in 2:nx-1
        shock_strength = abs(W.density_centers[i+1] - 2*W.density_centers[i] + W.density_centers[i-1])  # curvature of density
        flattening_coef = clamp(1.0 - 10.0 * shock_strength, 0.0, 1.0)
        slope_ρ[i] = slope_ρ[i] * flattening_coef
    end
end

function compute_flattening_coefficient(var::AbstractVector; ε=1e-10)
    n = length(var)
    θ = ones(n)  # start with no flattening (1.0)
    for i in 2:(n-1)
        Δp_plus = abs(var[i+1] - var[i])
        Δp_minus = abs(var[i] - var[i-1])
        Δp = max(Δp_plus, Δp_minus)
        p_avg = (abs(var[i+1]) + 2*abs(var[i]) + abs(var[i-1])) / 4

        # Avoid division by zero, ε small number
        shock_strength = Δp / (p_avg + ε)

        # Thresholds for flattening (tweak these!)
        if shock_strength > 0.1
            θ[i] = 0.0  # full flattening near strong shock
        elseif shock_strength > 0.05
            θ[i] = 0.5  # partial flattening near moderate gradient
        else
            θ[i] = 1.0  # no flattening in smooth regions
        end
    end

    # Boundaries: no flattening (or replicate neighbors)
    θ[1] = θ[2]
    θ[end] = θ[end-1]

    return θ
end


# This method finds a region, it needs to be tested
function shock_flattener(ρ, u, p, i, dx; tau=0.33, C=0.5)
    # First derivative of velocity (divergence/compression)
    divu = (u[i+1] - u[i-1]) / (2dx)

    # Pressure jump across stencil
    dp = p[i+1] - p[i-1]
    pmin = min(p[i-1], p[i+1], p[i]) + 1e-12  # avoid /0

    # Detect strong shock
    strong_shock = abs(dp) / pmin > tau
    compressive = divu < 0.0

    φ = 1.0
    if strong_shock && compressive
        # Scale coefficient: more flattening for stronger shocks
        φ = clamp(1.0 - C * abs(dp) / (p[i+1] + p[i-1] + 1e-12), 0.0, 1.0)
    end
    return φ
end
