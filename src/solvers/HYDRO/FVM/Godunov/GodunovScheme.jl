




function GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    

    Godunov_Log = nothing

    U_old = deepcopy(U)
    if :Debug ∉ features
        if :Verbose ∈ features
            Godunov_Log = open("Godunov_Log.txt", "a")
            write_solver_input(Godunov_Log, W, U_old, F)
        end

        if reconstruction != :Constant

            # Reconstruct left and right states at interfaces
            if reconstruction == :Linear

                slope_ρ = compute_slopes(W.density_centers, limiter)
                slope_u = compute_slopes(W.velocity_centers, limiter)
                slope_p = compute_slopes(W.pressure_centers, limiter)
                ρL, ρR = linear_reconstruct(W.density_centers, slope_ρ)
                uL, uR = linear_reconstruct(W.velocity_centers, slope_u)
                pL, pR = linear_reconstruct(W.pressure_centers, slope_p)

            elseif reconstruction == :Parabolic

                ρL, ρR = parabolic_reconstruct(W.density_centers)
                uL, uR = parabolic_reconstruct(W.velocity_centers)
                pL, pR = parabolic_reconstruct(W.pressure_centers)
            
            elseif reconstruction == :Cubic

                ρL, ρR = reconstruct_interfaces(W.density_centers, grid_points)
                uL, uR = reconstruct_interfaces(W.velocity_centers, grid_points)
                pL, pR = reconstruct_interfaces(W.pressure_centers, grid_points)

            end

            if flattening == true
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, ρL, ρR)
                #flat_coeff = compute_flattening_coefficient(W.pressure_centers)
                ρL, ρR = W.density_centers .+ flat_coeff .* (ρL .- W.density_centers), W.density_centers .+ flat_coeff .* (ρR .- W.density_centers)
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, uL, uR)
                uL, uR = W.velocity_centers .+ flat_coeff .* (uL .- W.velocity_centers), W.velocity_centers .+ flat_coeff .* (uR .- W.velocity_centers)
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, pL, pR)
                pL, pR = W.pressure_centers .+ flat_coeff .* (pL .- W.pressure_centers), W.pressure_centers .+ flat_coeff .* (pR .- W.pressure_centers)
            end
            
            if steepening == true
                steep_coeff = compute_steepening_coefficient(W.density_centers, W.pressure_centers)
                for i in 2:total_zones
                    δρ=ρR[i] - ρL[i]
                    ρL[i] = W.density_centers[i] - ((1+steep_coeff[i])/2) * δρ
                    ρR[i] = W.density_centers[i] + ((1+steep_coeff[i])/2) * δρ
                end
                #ρL, ρR = steepening_contact!(W.density_centers, W.pressure_centers, grid_points, γ, ρL, ρR)
            end

            #if reconstruction == :Parabolic || reconstruction == :Cubic 
            #    ρL, ρR = enforce_monotonicity!(W.density_centers, ρL, ρR)
            #    uL, uR = enforce_monotonicity!(W.velocity_centers, uL, uR)
            #    pL, pR = enforce_monotonicity!(W.pressure_centers, pL, pR)
            #end
        else
            ρL, ρR = W.density_centers, W.density_centers
            uL, uR = W.velocity_centers, W.velocity_centers
            pL, pR = W.pressure_centers, W.pressure_centers
        end

        apply_boundary_conditions(boundary_condition, U, zones, ghost_zones)
        
        if mode == :Standard

            for i in 2:total_zones
                @inbounds begin
                    if riemanntype == :HLL
                        f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                    elseif riemanntype == :HLLC
                        f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                    elseif riemanntype == :Exact
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

            for i in 2:total_zones-1
                @inbounds begin
                    U.density_centers[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                    U.momentum_centers[i] -= dt/spacing * (F.momentum_flux[i+1] - F.momentum_flux[i])
                    U.total_energy_centers[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])

                end
            end

        elseif mode == :Parallel
            Threads.@threads for i in 2:total_zones
                @inbounds begin
                    if riemanntype == :HLL
                        f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                    elseif riemanntype == :HLLC
                        f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                    elseif riemanntype == :Exact
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
            Threads.@threads for i in 2:total_zones-1
                @inbounds begin
                    U.density_centers[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                    U.momentum_centers[i] -= dt/spacing * (F.momentum_flux[i+1] - F.momentum_flux[i])
                    U.total_energy_centers[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])

                end
            end

        elseif mode == :GPU

            @warn "GPU mode not yet implemented for Godunov Schemes... Release coming in v0.2... Defaulting to Parallel..."
            GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for Lax Friedrichs... Release coming in v0.5... Defaulting to Parallel..."
            GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    
        end    

        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, cfl, spacing, dt; logfile=Godunov_Log)
            write_solver_output(Godunov_Log, W, U, F)
            close(Godunov_Log)
        end
    elseif :Debug ∈ features
        GodunovStep_Debug!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones)    
    end

end








function GodunovStepIx!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    

    Godunov_Log = nothing

    U_old = deepcopy(U)
    if :Debug ∉ features
        if :Verbose ∈ features
            Godunov_Log = open("Godunov_Log.txt", "a")
            write_solver_input(Godunov_Log, W, U_old, F)
        end

        if reconstruction != :Constant

            # Reconstruct left and right states at interfaces
            if reconstruction == :Linear

                slope_ρ = compute_slopes(W.density_centers, limiter)
                slope_u = compute_slopes(W.velocity_centers, limiter)
                slope_p = compute_slopes(W.pressure_centers, limiter)
                ρL, ρR = linear_reconstruct(W.density_centers, slope_ρ)
                uL, uR = linear_reconstruct(W.velocity_centers, slope_u)
                pL, pR = linear_reconstruct(W.pressure_centers, slope_p)

            elseif reconstruction == :Parabolic

                ρL, ρR = parabolic_reconstruct(W.density_centers)
                uL, uR = parabolic_reconstruct(W.velocity_centers)
                pL, pR = parabolic_reconstruct(W.pressure_centers)
            
            elseif reconstruction == :Cubic

                ρL, ρR = reconstruct_interfaces(W.density_centers, grid_points)
                uL, uR = reconstruct_interfaces(W.velocity_centers, grid_points)
                pL, pR = reconstruct_interfaces(W.pressure_centers, grid_points)

            end

            if flattening == true
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, ρL, ρR)
                #flat_coeff = compute_flattening_coefficient(W.pressure_centers)
                ρL, ρR = W.density_centers .+ flat_coeff .* (ρL .- W.density_centers), W.density_centers .+ flat_coeff .* (ρR .- W.density_centers)
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, uL, uR)
                uL, uR = W.velocity_centers .+ flat_coeff .* (uL .- W.velocity_centers), W.velocity_centers .+ flat_coeff .* (uR .- W.velocity_centers)
                flat_coeff = flatten_shocks!(W.density_centers, W.pressure_centers, W.velocity_centers, pL, pR)
                pL, pR = W.pressure_centers .+ flat_coeff .* (pL .- W.pressure_centers), W.pressure_centers .+ flat_coeff .* (pR .- W.pressure_centers)
            end
            
            if steepening == true
                #steep_coeff = compute_steepening_coefficient(W.density_centers, W.pressure_centers)
                #for i in 2:total_zones
                #    δρ=ρR[i] - ρL[i]
                #    ρL[i] = W.density_centers[i] - ((1+steep_coeff[i])/2) * δρ
                #    ρR[i] = W.density_centers[i] + ((1+steep_coeff[i])/2) * δρ
                #end
                ρL, ρR = steepening_contact!(W.density_centers, W.pressure_centers, grid_points, γ, ρL, ρR)
            end

            #if reconstruction == :Parabolic || reconstruction == :Cubic 
            #    ρL, ρR = enforce_monotonicity!(W.density_centers, ρL, ρR)
            #    uL, uR = enforce_monotonicity!(W.velocity_centers, uL, uR)
            #    pL, pR = enforce_monotonicity!(W.pressure_centers, pL, pR)
            #end
        else
            ρL, ρR = W.density_centers, W.density_centers
            uL, uR = W.velocity_centers, W.velocity_centers
            pL, pR = W.pressure_centers, W.pressure_centers
        end

        apply_boundary_conditions(boundary_condition, U, zones, ghost_zones)
        
        if mode == :Standard

            for i in 2:total_zones
                @inbounds begin
                    if riemanntype == :HLL
                        f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                    elseif riemanntype == :HLLC
                        f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                    elseif riemanntype == :Exact
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

            for i in 2:total_zones-1
                @inbounds begin
                    U.density_centers[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                    U.momentum_centers[i] -= dt/spacing * (F.momentum_flux[i+1] - F.momentum_flux[i])
                    U.total_energy_centers[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])

                end
            end

        elseif mode == :Parallel
            Threads.@threads for i in 2:total_zones
                @inbounds begin
                    if riemanntype == :HLL
                        f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)
                    elseif riemanntype == :HLLC
                        f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
                    elseif riemanntype == :Exact
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
            Threads.@threads for i in 2:total_zones-1
                @inbounds begin
                    U.density_centers[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                    U.momentum_centers[i] -= dt/spacing * (F.momentum_flux[i+1] - F.momentum_flux[i])
                    U.total_energy_centers[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])

                end
            end

        elseif mode == :GPU

            @warn "GPU mode not yet implemented for Godunov Schemes... Release coming in v0.2... Defaulting to Parallel..."
            GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for Lax Friedrichs... Release coming in v0.5... Defaulting to Parallel..."
            GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones, grid_points)    
        end    

        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, cfl, spacing, dt; logfile=Godunov_Log)
            write_solver_output(Godunov_Log, W, U, F)
            close(Godunov_Log)
        end
    elseif :Debug ∈ features
        GodunovStep_Debug!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, total_zones, zones, ghost_zones)    
    end

end