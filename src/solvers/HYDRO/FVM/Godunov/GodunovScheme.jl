function GodunovStep!(
    W, U, F,
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

    # --- Flux computation ---
    if mode == :Standard
        for i in 2:total_zones
            @inbounds begin
                f1 = f2 = f3 = 0.0
                if riemanntype == :HLL
                    f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                elseif riemanntype == :HLLC
                    f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                elseif riemanntype == :Exact
                    pstar, ustar = solve_star_pressure(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
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

        for i in 2:total_zones-1
            @inbounds begin
                U.centers.density[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                U.centers.momentum[1][i] -= dt/spacing * (F.momentum_flux[1][i+1] - F.momentum_flux[1][i])
                U.centers.total_energy[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])
            end
        end

    elseif mode == :Parallel
        Threads.@threads for i in 2:total_zones
            @inbounds begin
                f1 = f2 = f3 = 0.0
                if riemanntype == :HLL
                    f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                elseif riemanntype == :HLLC
                    f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                elseif riemanntype == :Exact
                    pstar, ustar = solve_star_pressure(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
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

        Threads.@threads for i in 2:total_zones-1
            @inbounds begin
                U.centers.density[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                U.centers.momentum[1][i] -= dt/spacing * (F.momentum_flux[1][i+1] - F.momentum_flux[1][i])
                U.centers.total_energy[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])
            end
        end
    end
end
