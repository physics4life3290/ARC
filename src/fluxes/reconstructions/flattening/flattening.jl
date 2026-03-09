


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

# ============================================================================= #
#                                                                               #
#   Colella-Woodward Flattening Coefficient                                    #
#   Algorithms for Realistic Computations (ARC)                                #
#   Reference: Colella & Woodward (1984), JCP 54, 174                          #
#                                                                               #
# ============================================================================= #

"""
    compute_flattening_coeff!(Fi, P, u; z1, z2, δp_threshold, ϵ)

Compute the Colella-Woodward flattening coefficient in-place.

Fi[i] ∈ [0,1] where:
  0 → no flattening (full high-order reconstruction)
  1 → full flattening (drop to piecewise constant)

Arguments:
  Fi            — output flattening coefficient vector (modified in-place)
  P             — pressure array
  u             — velocity array
  z1            — lower transition zone boundary  (default 0.75, C-W standard)
  z2            — upper transition zone boundary  (default 0.85, C-W standard)
  δp_threshold  — pressure jump threshold         (default 0.33, C-W standard)
  ϵ             — small number to prevent division by zero
  widen         — number of widening passes       (default 2)
"""
function compute_flattening_coeff!(
    Fi::Vector{Float64},
    P::Vector{Float64},
    u::Vector{Float64},
    Fi_tmp::Vector{Float64};          # pass pre-allocated buffer — no allocation inside
    z1::Float64            = 0.75,
    z2::Float64            = 0.85,
    δp_threshold::Float64  = 0.33,
    ϵ::Float64             = 1.0e-12,
    widen::Int             = 2
)
    n = length(P)
    fill!(Fi, 0.0)

    # ----------------------------------------------------------------- #
    #   Step 1: Compute raw flattening at each cell                      #
    # ----------------------------------------------------------------- #
    @inbounds for i in 3:n-2
        num   = abs(P[i+1] - P[i-1])
        den   = abs(P[i+2] - P[i-2]) + ϵ
        ratio = num / den

        # Colella-Woodward transition zone [z1, z2]
        fbar = if ratio <= z1
            0.0
        elseif ratio >= z2
            1.0
        else
            (ratio - z1) / (z2 - z1)
        end

        # Only flatten at shocks: strong pressure jump AND compressive flow
        Pmin          = min(P[i+1], P[i-1])
        pressure_jump = (num / (abs(Pmin) + ϵ)) > δp_threshold
        compressive   = (u[i+1] - u[i-1]) < 0.0

        Fi[i] = (pressure_jump && compressive) ? fbar : 0.0
    end

    # ----------------------------------------------------------------- #
    #   Step 2: Widen the flattening zone                                #
    #   Multiple passes ensure the shock peak and its neighbors are      #
    #   all flattened — prevents odd-even decoupling at the shock front  #
    # ----------------------------------------------------------------- #
    for _ in 1:widen
        copyto!(Fi_tmp, Fi)
        @inbounds for i in 3:n-2
            # Spread toward the high-pressure side
            Fi[i] = if P[i+1] < P[i-1]
                max(Fi_tmp[i], Fi_tmp[i+1])
            else
                max(Fi_tmp[i], Fi_tmp[i-1])
            end
        end
    end

    return Fi
end


# ============================================================================= #
#   APPLICATION: Flatten interface states                                        #
#                                                                               #
#   Given left/right reconstructed states at interface i+1/2,                  #
#   blend toward piecewise constant using flattening coefficient Fi.            #
#                                                                               #
#   Q_L_flat = (1 - Fi) * Q_L + Fi * Q[i]                                     #
#   Q_R_flat = (1 - Fi) * Q_R + Fi * Q[i+1]                                   #
# ============================================================================= #

@inline function apply_flattening(QL::Float64, QR::Float64,
                                   Qi::Float64,  Qi1::Float64,
                                   Fi::Float64)
    QL_flat = (1.0 - Fi) * QL + Fi * Qi
    QR_flat = (1.0 - Fi) * QR + Fi * Qi1
    return QL_flat, QR_flat
end