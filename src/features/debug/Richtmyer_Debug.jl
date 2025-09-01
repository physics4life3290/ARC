




function RichtmyerStep_Debug!(W::PrimitiveVariables, U::ConservativeVariables,
                        F::FluxVariables,
                        _grid::CartesianGrid1D,
                        UserInput,
                        dt::Float64,
                        ghost_zones, total_zones, spacing, γ, mode, features, cfl)
    
    try
        U_old = deepcopy(U)
        CalculateFlux!(W, U,F)
        if :Verbose ∈ features
                
            Richtmyer_Log = open("Richtmyer_Log.txt", "a")
            write_solver_input(Richtmyer_Log, W, U_old, F)
        end
        # Allocate U_half as a separate copy to avoid overwriting U
        U_half = ConservativeVariables(
            copy(U.density_centers),
            copy(U.momentum_centers),
            copy(U.total_energy_centers),
            nothing,
            nothing, 
            nothing
        )

        # Start modes here #
        if mode == :Standard

            for i in (ghost_zones + 1):(total_zones - ghost_zones)
                @inbounds begin
                    ic, il = i, i-1 

                    U_half.density_centers[i] = 0.5*(U.density_centers[ic] + U.density_centers[il]) -
                                                (dt/(2*spacing))*(F.density_flux[ic] - F.density_flux[il])
                    U_half.momentum_centers[i] = 0.5*(U.momentum_centers[ic] + U.momentum_centers[il]) -
                                                (dt/(2*spacing))*(F.momentum_flux[ic] - F.momentum_flux[il])
                    U_half.total_energy_centers[i] = 0.5*(U.total_energy_centers[ic] + U.total_energy_centers[il]) -
                                                    (dt/(2*spacing))*(F.total_energy_flux[ic] - F.total_energy_flux[il])
                end
            end

        elseif mode == :Parallel
            # Predictor step
            Threads.@threads for i in (ghost_zones + 1):(total_zones - ghost_zones)
                @inbounds begin
                    ic, il = i, i-1 

                    U_half.density_centers[i] = 0.5*(U.density_centers[ic] + U.density_centers[il]) -
                                                (dt/(2*spacing))*(F.density_flux[ic] - F.density_flux[il])
                    U_half.momentum_centers[i] = 0.5*(U.momentum_centers[ic] + U.momentum_centers[il]) -
                                                (dt/(2*spacing))*(F.momentum_flux[ic] - F.momentum_flux[il])
                    U_half.total_energy_centers[i] = 0.5*(U.total_energy_centers[ic] + U.total_energy_centers[il]) -
                                                    (dt/(2*spacing))*(F.total_energy_flux[ic] - F.total_energy_flux[il])
                end
            end
        elseif mode == :GPU

            @warn "GPU mode not yet implemented for Richtmyer Method... Release coming in v0.2... Defaulting to Parallel..."
            Richtmyer_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for Richtmyer Method... Release coming in v0.5... Defaulting to Parallel..."
            Richtmyer_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        end

        # Reconstruct primitive variables from U_half
        W_half = PrimitiveVariables(
            copy(U_half.density_centers),
            U_half.momentum_centers ./ U_half.density_centers,
            zeros(length(U_half.density_centers)),
            zeros(length(U_half.density_centers)),
            nothing, nothing, nothing, nothing
        )
        W_half.pressure_centers .= (γ - 1) .* (U_half.total_energy_centers .- 0.5 .* U_half.density_centers .* W_half.velocity_centers .^ 2)
        W_half.internal_energy_centers .= W_half.pressure_centers ./ ((γ-1) .* U_half.density_centers)

        apply_boundary_conditions(UserInput.Primary_Input.boundary_condition, U_half, _grid.coord1.zones, _grid.coord1.ghost_zones)

        # Compute F_half from U_half
        F_half = FluxVariables(
            U_half.momentum_centers,
            U_half.momentum_centers .* W_half.velocity_centers .+ W_half.pressure_centers,
            W_half.velocity_centers .* (U_half.total_energy_centers .+ W_half.pressure_centers)
        )
    
        if mode == :Standard

            for i in (ghost_zones+1):(total_zones-ghost_zones)
                @inbounds begin
                    U.density_centers[i]      -= (dt/spacing)*(F_half.density_flux[i+1] - F_half.density_flux[i])
                    U.momentum_centers[i]     -= (dt/spacing)*(F_half.momentum_flux[i+1] - F_half.momentum_flux[i])
                    U.total_energy_centers[i] -= (dt/spacing)*(F_half.total_energy_flux[i+1] - F_half.total_energy_flux[i])
                end
            end

        elseif mode == :Parallel

            # Corrector step
            Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
                @inbounds begin
                    U.density_centers[i]      -= (dt/spacing)*(F_half.density_flux[i+1] - F_half.density_flux[i])
                    U.momentum_centers[i]     -= (dt/spacing)*(F_half.momentum_flux[i+1] - F_half.momentum_flux[i])
                    U.total_energy_centers[i] -= (dt/spacing)*(F_half.total_energy_flux[i+1] - F_half.total_energy_flux[i])
                end
            end

        elseif mode == :GPU

            @warn "GPU mode not yet implemented for Lax Friedrichs... Release coming in v0.2... Defaulting to Parallel..."
            RichtmyerStep!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for Lax Friedrichs... Release coming in v0.5... Defaulting to Parallel..."
            RichtmyerStep!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        end

        # Final boundary conditions on updated U
        apply_boundary_conditions(UserInput.Primary_Input.boundary_condition, U, _grid.coord1.zones, _grid.coord1.ghost_zones)
        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, cfl, spacing, dt; logfile=Richtmyer_Log)
            write_solver_output(Richtmyer_Log, W, U, F)
            close(Richtmyer_Log)
        end 
    catch e 
        println("The incoming data is 
        Density: $(U_old.density_centers)
        Length of Density: $(length(U_old.density_centers))
        Momentum: $(U_old.momentum_centers)
        Length of Momentum: $(length(U_old.momentum_centers))
        Total Energy: $(U_old.total_energy_centers)
        Length of Total Energy: $(length(U_old.total_energy_centers))")

        @error "An error occured during Richtmyer step: $e"
    end
end