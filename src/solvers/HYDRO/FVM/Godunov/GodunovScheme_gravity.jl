
function GodunovStep!(
    W, U, F, S,
    reconstruction, limiter,
    flattening, steepening,
    boundary_condition, riemanntype,
    γ, spacing, dt,
    total_zones, zones, ghost_zones, grid_points,
    mode
)

    # --- Reconstruction ---
    if reconstruction != :Constant
        if reconstruction == :Linear
            slope_ρ = compute_slopes(W.centers.density, limiter)
            slope_u = compute_slopes(W.centers.velocity[1], limiter)
            slope_p = compute_slopes(W.centers.pressure, limiter)
            ρL, ρR = linear_reconstruct(W.centers.density, slope_ρ)
            uL, uR = linear_reconstruct(W.centers.velocity[1], slope_u)
            pL, pR = linear_reconstruct(W.centers.pressure, slope_p)

        elseif reconstruction == :Parabolic
            ρL, ρR = parabolic_reconstruct(W.centers.density)
            uL, uR = parabolic_reconstruct(W.centers.velocity[1])
            pL, pR = parabolic_reconstruct(W.centers.pressure)

        elseif reconstruction == :Cubic
            ρL, ρR = reconstruct_interfaces(W.centers.density, grid_points)
            uL, uR = reconstruct_interfaces(W.centers.velocity[1], grid_points)
            pL, pR = reconstruct_interfaces(W.centers.pressure, grid_points)
        end

        if steepening
            ρL, ρR = steepening_contact!(W.centers.density, W.centers.pressure, grid_points, γ, ρL, ρR)
        end

        if flattening
            flat_coeff = compute_flattening_coeff!(zeros(length(W.centers.density)), W.centers.pressure, W.centers.velocity[1])
            if reconstruction == :Linear
                @inbounds begin
                    ρL[1:end-1] .= W.centers.density .+ flat_coeff .* (ρL[1:end-1] .- W.centers.density)
                    ρR[1:end-1] .= W.centers.density .+ flat_coeff .* (ρR[1:end-1] .- W.centers.density)
                    uL[1:end-1] .= W.centers.velocity[1] .+ flat_coeff .* (uL[1:end-1] .- W.centers.velocity[1])
                    uR[1:end-1] .= W.centers.velocity[1] .+ flat_coeff .* (uR[1:end-1] .- W.centers.velocity[1])
                    pL[1:end-1] .= W.centers.pressure .+ flat_coeff .* (pL[1:end-1] .- W.centers.pressure)
                    pR[1:end-1] .= W.centers.pressure .+ flat_coeff .* (pR[1:end-1] .- W.centers.pressure)
                end
            else
                @inbounds begin
                    ρL .= W.centers.density .+ flat_coeff .* (ρL .- W.centers.density)
                    ρR .= W.centers.density .+ flat_coeff .* (ρR .- W.centers.density)
                    uL .= W.centers.velocity[1] .+ flat_coeff .* (uL .- W.centers.velocity[1])
                    uR .= W.centers.velocity[1] .+ flat_coeff .* (uR .- W.centers.velocity[1])
                    pL .= W.centers.pressure .+ flat_coeff .* (pL .- W.centers.pressure)
                    pR .= W.centers.pressure .+ flat_coeff .* (pR .- W.centers.pressure)
                end
            end
    
        end
    else
        ρL, ρR = W.centers.density, W.centers.density
        uL, uR = W.centers.velocity[1], W.centers.velocity[1]
        pL, pR = W.centers.pressure, W.centers.pressure
    end

    # --- Boundary conditions ---
    apply_boundary_conditions(boundary_condition, U, zones, ghost_zones)

    # --- 1. Compute Gravitational Acceleration ---
    # g = -dΦ/dr
    g = zeros(eltype(S), total_zones)
    for i in 2:total_zones-1
        # Central difference for acceleration at cell centers
        g[i] = -(S[i+1] - S[i-1]) / (2.0 * spacing)
    end
    # Simple extrapolation for boundaries
    g[1] = -(S[2] - S[1]) / spacing
    g[total_zones] = -(S[total_zones] - S[total_zones-1]) / spacing
    p0 = -W.centers.density .* g .* spacing
    # --- 2. Flux computation ---
    if mode == :Standard || mode == :Parallel
        # (Select loop type based on mode as in your original code)
        # Calculate Fluxes F.density_flux, F.momentum_flux, F.total_energy_flux...
        for i in 2:total_zones
            @inbounds begin
                f1 = f2 = f3 = 0.0
                if riemanntype == :HLL
                    f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                elseif riemanntype == :HLLC
                    f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                elseif riemanntype == :Exact
                    pstar, ustar = solve_star_pressure(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ, p0=p0[i])
                    prims = sample_exact(0.0, ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], pstar, ustar, γ)
                    f1 = prims[1] * prims[2]
                    f2 = prims[1] * prims[2]^2 + prims[3]
                    E = prims[3]/(γ-1.0) + 0.5*prims[1]*prims[2]^2
                    f3 = prims[2] * (E + prims[3])
                end
                F.density_flux[i] = f1
                F.momentum_flux[1][i] = f2
                F.total_energy_flux[i] = f3
            end
        end

        # This part updates U using only the fluxes (Standard Conservation Law)
        for i in (ghost_zones + 1):(total_zones - ghost_zones)
            @inbounds begin
                dF1 = (F.density_flux[i+1] - F.density_flux[i]) / spacing
                dF2 = (F.momentum_flux[1][i+1] - F.momentum_flux[1][i]) / spacing
                dF3 = (F.total_energy_flux[i+1] - F.total_energy_flux[i]) / spacing

                U.centers.density[i]      -= dt * dF1
                U.centers.momentum[1][i]  -= dt * dF2
                U.centers.total_energy[i] -= dt * dF3
            end
        end

        # --- Step 2: Source Terms (Lie Splitting) ---
        # Update the conserved variables based on geometry and gravity
        for i in (ghost_zones + 1):(total_zones - ghost_zones)
            @inbounds begin
                # Get primitive values from the *updated* or *current* state 
                # (For Lie splitting, we use the intermediate U to get current prims)
                ri = grid_u.axes[1].centers[i]
                ρ  = U.centers.density[i]
                m  = U.centers.momentum[1][i]
                E  = U.centers.total_energy[i]
                
                # Derived primitives
                u  = m / ρ
                # p = (γ - 1) * (E - 0.5 * ρ * u^2)
                p  = (γ - 1.0) * (E - 0.5 * m^2 / ρ)

                # 1. Spherical Geometric terms (-2 * flux / r)
                # 2. Gravitational terms
                s_ρ = - (2.0 * ρ * u) / ri
                s_m = - (2.0 * ρ * u^2) / ri + ρ * g[i]
                s_E = - (2.0 * u * (E + p)) / ri + ρ * u * g[i]

                # Final Update
                U.centers.density[i]      += dt * s_ρ
                U.centers.momentum[1][i]  += dt * s_m
                U.centers.total_energy[i] += dt * s_E
            end
        end

    end
end