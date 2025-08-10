




# Minmod limiter
function minmod(a, b)
    if a*b <= 0
        return 0.0
    else
        return sign(a) * min(abs(a), abs(b))
    end
end

# I need to update all of this to match my ARC format instead of my prototype format.
# I also need to implement the parabolic and cubic interpolations via lagrange interpolation (This will prepare for a non uniform grid)

function ppm_reconstruction(U, _grid)
    nx = size(U, 2)
    UL = copy(U)
    UR = copy(U)

    for i in 2:nx-1
        for v in 1:3
            dL = U[v, i] - U[v, i-1]
            dR = U[v, i+1] - U[v, i]
            dC = 0.5 * (U[v, i+1] - U[v, i-1])

            # Minmod slope as base
            slope = minmod(dL, dR)

            # --------------------------
            # Flattening near strong shocks
            # --------------------------
            shock_strength = abs(U[1, i+1] - 2U[1, i] + U[1, i-1])  # curvature of density
            flattening_coef = clamp(1.0 - 10.0 * shock_strength, 0.0, 1.0)
            slope *= flattening_coef

            # --------------------------
            # Steepening for contact discontinuities
            # --------------------------
            if v == 1  # only steepen density
                ΔpL = abs(U[3, i] - U[3, i-1])
                ΔpR = abs(U[3, i+1] - U[3, i])
                ΔρL = abs(U[1, i] - U[1, i-1])
                ΔρR = abs(U[1, i+1] - U[1, i])
                steepening_coef = 0.0

                if (ΔρL > 0.02 * U[1, i]) && (ΔpL < 0.1 * U[3, i]) &&
                   (ΔρR > 0.02 * U[1, i]) && (ΔpR < 0.1 * U[3, i])
                    steepening_coef = 0.3  # strengthen slope slightly
                end

                slope += steepening_coef * dC
            end

            # Final interface values
            UL[v, i] = U[v, i] - 0.5 * slope
            UR[v, i] = U[v, i] + 0.5 * slope

            # Monotonicity constraint: Clamp to local extrema
            umin = min(U[v, i-1], U[v, i], U[v, i+1])
            umax = max(U[v, i-1], U[v, i], U[v, i+1])
            UL[v, i] = clamp(UL[v, i], umin, umax)
            UR[v, i] = clamp(UR[v, i], umin, umax)
        end
    end

    return UL, UR
end

function PPM_Step!(UserInput, _grid, W, U, t)
    nx = _grid.xcoord.total_zones
    γ = UserInput.secondary_input.γ
    dx = _grid.xcoord.spacing
    ρR, uR, pR = zeros(nx), zeros(nx), zeros(nx)
    ρL, uL, pL = zeros(nx), zeros(nx), zeros(nx)
    # Compute slopes for primitive variables
    
    apply_boundary_conditions(UserInput, U, _grid)
    # Reconstruct left and right states at interfaces
    UL, UR = ppm_reconstruction(U)

    for i in 1:nx
        ρR[i] = UR[1, i]
        uR[i] = UR[2, i] / UR[1, i]
        pR[i] = (γ - 1) * (UR[3, i] - 0.5 * UR[2, i]^2 / UR[1, i])

        ρL[i] = UL[1, i]
        uL[i] = UL[2, i] / UL[1, i]
        pL[i] = (γ - 1) * (UL[3, i] - 0.5 * UL[2, i]^2 / UL[1, i])
    end

    
    # Fluxes at interfaces
    F1 = zeros(nx+1)
    F2 = zeros(nx+1)
    F3 = zeros(nx+1)
    Threads.@threads for i in 2:nx
        @inbounds begin
            #f1, f2, f3 = riemann_exact(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i])
            f1, f2, f3 = ExactRiemannSolve!((ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i]),γ)
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