




function LaxFriedrichs_Step_Debug!(W::PrimitiveVariables, U::ConservativeVariables, F::FluxVariables, dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64, mode::Symbol, features::Vector{Symbol}, cfl)
    
    try
        LaxFriedrichs_Log = nothing
        # Store old state for increment norms
        U_old = deepcopy(U)
        # Calculate fluxes
        CalculateFlux!(W, U, F)

        if :Verbose ∈ features
            
            LaxFriedrichs_Log = open("LaxFriedrichs_Log.txt", "a")
            write_solver_input(LaxFriedrichs_Log, W, U_old, F)
        end
        
        if mode == :Standard

            for i in ghost_zones+1:total_zones-ghost_zones
                @inbounds begin
                    U.density_centers[i]      = 0.5*(U.density_centers[i+1] + U.density_centers[i-1]) -
                                                (dt/(2*spacing))*(F.density_flux[i+1]   - F.density_flux[i-1])
                    U.momentum_centers[i]     = 0.5*(U.momentum_centers[i+1] + U.momentum_centers[i-1]) -
                                                (dt/(2*spacing))*(F.momentum_flux[i+1]   - F.momentum_flux[i-1])
                    U.total_energy_centers[i] = 0.5*(U.total_energy_centers[i+1] + U.total_energy_centers[i-1]) -
                                                (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
                end
            end

        elseif mode == :Parallel
        
            Threads.@threads for i in ghost_zones+1:total_zones-ghost_zones
                @inbounds begin
                    U.density_centers[i]      = 0.5*(U.density_centers[i+1] + U.density_centers[i-1]) -
                                                (dt/(2*spacing))*(F.density_flux[i+1]   - F.density_flux[i-1])
                    U.momentum_centers[i]     = 0.5*(U.momentum_centers[i+1] + U.momentum_centers[i-1]) -
                                                (dt/(2*spacing))*(F.momentum_flux[i+1]   - F.momentum_flux[i-1])
                    U.total_energy_centers[i] = 0.5*(U.total_energy_centers[i+1] + U.total_energy_centers[i-1]) -
                                                (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
                end
            end

        elseif mode == :GPU

            @warn "GPU mode not yet implemented for Lax Friedrichs... Release coming in v0.2... Defaulting to Parallel..."
            LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for Lax Friedrichs... Release coming in v0.5... Defaulting to Parallel..."
            LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        end

        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, cfl, spacing, dt; logfile=LaxFriedrichs_Log)
            write_solver_output(LaxFriedrichs_Log, W, U, F)
            close(LaxFriedrichs_Log)
        end 
    catch e 
        println("The incoming data is 
        Density: $(U_old.density_centers)
        Length of Density: $(length(U_old.density_centers))
        Momentum: $(U_old.momentum_centers)
        Length of Momentum: $(length(U_old.momentum_centers))
        Total Energy: $(U_old.total_energy_centers)
        Length of Total Energy: $(length(U_old.total_energy_centers))")

        @error "An error occured during Lax Friedrichs step: $e"
    end

end