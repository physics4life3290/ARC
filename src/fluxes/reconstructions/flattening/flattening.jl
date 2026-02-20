


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

function compute_flattening_coefficient_2d(var::AbstractMatrix; ε=1e-10)
    nx, ny = size(var)
    θ = ones(nx, ny)  # start with no flattening (1.0)

    for j in 2:(ny-1)
        for i in 2:(nx-1)
            @inbounds begin
                # --- X-direction ---
                Δx_plus  = abs(var[i+1,j] - var[i,j])
                Δx_minus = abs(var[i,j] - var[i-1,j])
                Δx = max(Δx_plus, Δx_minus)
                p_avg_x = (abs(var[i+1,j]) + 2*abs(var[i,j]) + abs(var[i-1,j])) / 4
                shock_x = Δx / (p_avg_x + ε)

                # --- Y-direction ---
                Δy_plus  = abs(var[i,j+1] - var[i,j])
                Δy_minus = abs(var[i,j] - var[i,j-1])
                Δy = max(Δy_plus, Δy_minus)
                p_avg_y = (abs(var[i,j+1]) + 2*abs(var[i,j]) + abs(var[i,j-1])) / 4
                shock_y = Δy / (p_avg_y + ε)

                # Combine X and Y contributions (take max)
                shock_strength = max(shock_x, shock_y)

                # Flattening thresholds (tweakable)
                if shock_strength > 0.1
                    θ[i,j] = 0.0  # full flattening near strong shock
                elseif shock_strength > 0.05
                    θ[i,j] = 0.5  # partial flattening
                else
                    θ[i,j] = 1.0  # smooth region
                end
            end
        end
    end

    # Boundaries: replicate nearest interior cell
    θ[1,:] .= θ[2,:]
    θ[end,:] .= θ[end-1,:]
    θ[:,1] .= θ[:,2]
    θ[:,end] .= θ[:,end-1]

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

# --- Flattening coefficient ---
function compute_flattening_coeff!(
    Fi::Vector{Float64},
    P::Vector{Float64},
    u::Vector{Float64};
    ϵ = 1e-12
)
    n = length(P)
    fill!(Fi, 0.0)
    @inbounds for i in 3:n-2
        num = abs(P[i+1] - P[i-1])
        den = abs(P[i+2] - P[i-2]) + ϵ
        ratio = num / den
        fbar = clamp(4.0*(ratio-0.5),0.0,1.0)
        Pmin = max(min(P[i+1],P[i-1]),ϵ)
        pressure_jump = (num/Pmin) > 0.33
        compressive = (u[i+1]-u[i-1]) < 0.0
        Fi[i] = (pressure_jump && compressive) ? fbar : 0.0
    end

    # Widening
    Fi_tmp = copy(Fi)
    @inbounds for i in 3:n-2
        Fi[i] = P[i+1] < P[i-1] ? max(Fi_tmp[i], Fi_tmp[i+1]) : max(Fi_tmp[i], Fi_tmp[i-1])
    end

    return Fi
end