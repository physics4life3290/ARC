




function GodunovStep!(user_input, _grid, W, U, F, dt)
    
    nx = _grid.coord1.total_zones
    γ = user_input.Secondary_Input.gamma
    dx = _grid.coord1.spacing

    if user_input.Solver_Input.reconstruction != :Constant

        # Reconstruct left and right states at interfaces
        if user_input.Solver_Input.reconstruction == :Linear
            slope_ρ = compute_slopes(W.density_centers, user_input.Solver_Input.limiter)
            slope_u = compute_slopes(W.velocity_centers, user_input.Solver_Input.limiter)
            slope_p = compute_slopes(W.pressure_centers, user_input.Solver_Input.limiter)
            ρL, ρR = linear_reconstruct(W.density_centers, slope_ρ)
            uL, uR = linear_reconstruct(W.velocity_centers, slope_u)
            pL, pR = linear_reconstruct(W.pressure_centers, slope_p)
        elseif user_input.Solver_Input.reconstruction == :Parabolic
            ρL, ρR = parabolic_reconstruct(W.density_centers)
            uL, uR = parabolic_reconstruct(W.velocity_centers)
            pL, pR = parabolic_reconstruct(W.pressure_centers)
            #if user_input.Solver_Input.limiter !== nothing 
            #    ρL, ρR = apply_ppm_limiter!(W.density_centers, ρL, ρR, user_input)
            #    uL, uR = apply_ppm_limiter!(W.velocity_centers, uL, uR, user_input)
            #    pL, pR = apply_ppm_limiter!(W.pressure_centers, pL,pR, user_input)
            #end
        end

        if user_input.Solver_Input.flattening == true
            flat_coeff = compute_flattening_coefficient(W.pressure_centers)
            ρL, ρR = W.density_centers .+ flat_coeff .* (ρL .- W.density_centers), W.density_centers .+ flat_coeff .* (ρR .- W.density_centers)
            uL, uR = W.velocity_centers .+ flat_coeff .* (uL .- W.velocity_centers), W.velocity_centers .+ flat_coeff .* (uR .- W.velocity_centers)
            pL, pR = W.pressure_centers .+ flat_coeff .* (pL .- W.pressure_centers), W.pressure_centers .+ flat_coeff .* (pR .- W.pressure_centers)
        end
        
        if user_input.Solver_Input.steepening == true
            steep_coeff = compute_steepening_coefficient(W.density_centers, W.pressure_centers)
            for i in 2:nx
                δρ=ρR[i] - ρL[i]
                ρL[i] = W.density_centers[i] - ((1+steep_coeff[i])/2) * δρ
                ρR[i] = W.density_centers[i] + ((1+steep_coeff[i])/2) * δρ
            end
        
        end

    else
        ρL, ρR = W.density_centers, W.density_centers
        uL, uR = W.velocity_centers, W.velocity_centers
        pL, pR = W.pressure_centers, W.pressure_centers
    end

    apply_boundary_conditions(user_input, U, _grid)
    # Fluxes at interfaces
    F1 = zeros(nx+1)
    F2 = zeros(nx+1)
    F3 = zeros(nx+1)
    Threads.@threads for i in 2:nx
        @inbounds begin
            if user_input.Solver_Input.riemanntype == :HLL
                f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
            elseif user_input.Solver_Input.riemanntype == :HLLC
                f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
            elseif user_input.Solver_Input.riemanntype == :Exact
                #f1, f2, f3 = ExactRiemannSolve!((ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i]),γ)
                pstar, ustar = solve_star_pressure(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                prims = sample_exact(0.0, ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], pstar, ustar, γ)
                f1 = prims[1] * prims[2]
                f2 = prims[1] * prims[2]^2 + prims[3]
                E = prims[3] / (γ - 1.0) + 0.5 * prims[1] * prims[2]^2
                f3 = prims[2] * (E + prims[3])
            end
            F.density_flux[i] = f1
            F.momentum_flux[i] = f2
            F.total_energy_flux[i] = f3
        end
    end

    # Update conserved variables (interior cells)
    Threads.@threads for i in 2:nx-1
        @inbounds begin
            U.density_centers[i] -= dt/dx * (F.density_flux[i+1] - F.density_flux[i])
            U.momentum_centers[i] -= dt/dx * (F.momentum_flux[i+1] - F.momentum_flux[i])
            U.total_energy_centers[i] -= dt/dx * (F.total_energy_flux[i+1] - F.total_energy_flux[i])

        end
    end
end