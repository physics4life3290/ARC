

function LaxFriedrichs_Step!(W::PrimitiveVariables, U::ConservativeVariables, F::FluxVariables, dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64, mode::Symbol, features::Vector{Symbol}, cfl)
    
    LaxFriedrichs_Log = nothing
    # Store old state for increment norms
    U_old = deepcopy(U)
    # Calculate fluxes
    CalculateFlux!(W, U, F)
    if :Debug ∉ features
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
    elseif :Debug ∈ features
        LaxFriedrichs_Step_Debug!(W, U, F, dt, ghost_zones, total_zones, spacing, mode, features, cfl)
    end

end

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


#= From here down is still in testing =#

function LaxFriedrichs_ViscousStep!(W, U::ConservativeVariables, F::FluxVariables,
                                    dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64;
                                    alpha_visc::Float64 = 0.5,
                                    rhs_density = nothing, rhs_momentum = nothing, rhs_energy = nothing)

    # ensure fluxes are up-to-date
    CalculateFlux!(W, U, F)

    ncells = total_zones
    # allocate or validate RHS buffers (one-per-cell for the updated value)
    if rhs_density === nothing
        rhs_density = similar(U.density_centers)
    else
        @assert length(rhs_density) == length(U.density_centers)
    end
    if rhs_momentum === nothing
        rhs_momentum = similar(U.momentum_centers)
    else
        @assert length(rhs_momentum) == length(U.momentum_centers)
    end
    if rhs_energy === nothing
        rhs_energy = similar(U.total_energy_centers)
    else
        @assert length(rhs_energy) == length(U.total_energy_centers)
    end

    invdx = dt/spacing   # used for flux divergence multiplication

    # Compute updates into rhs_* (thread-safe: writes to unique i indices; reads only U and F)
    Threads.@threads for i in (ghost_zones+1):(total_zones - ghost_zones)
        @inbounds begin
            # flux divergence (central two-point form used in your original)
            # (0.5*(F[i]+F[i+1]) - 0.5*(F[i-1]+F[i])) = 0.5*(F[i+1]-F[i-1])
            fluxdiv_density  = 0.5*(F.density_flux[i+1]  - F.density_flux[i-1])
            fluxdiv_momentum = 0.5*(F.momentum_flux[i+1] - F.momentum_flux[i-1])
            fluxdiv_energy   = 0.5*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])

            # laplacian (second difference) for artificial viscosity
            lap_density  = (U.density_centers[i+1]      - 2.0*U.density_centers[i]      + U.density_centers[i-1])
            lap_momentum = (U.momentum_centers[i+1]     - 2.0*U.momentum_centers[i]     + U.momentum_centers[i-1])
            lap_energy   = (U.total_energy_centers[i+1] - 2.0*U.total_energy_centers[i] + U.total_energy_centers[i-1])

            # explicit update: U_new = U_old - (dt/dx)*fluxdiv + alpha_visc * lap
            rhs_density[i]  = U.density_centers[i]      - invdx*fluxdiv_density  + alpha_visc * lap_density
            rhs_momentum[i] = U.momentum_centers[i]     - invdx*fluxdiv_momentum + alpha_visc * lap_momentum
            rhs_energy[i]   = U.total_energy_centers[i] - invdx*fluxdiv_energy   + alpha_visc * lap_energy
        end
    end

    # copy updated values back into U for the active cells (thread-safe: distinct cell writes)
    Threads.@threads for i in (ghost_zones+1):(total_zones - ghost_zones)
        @inbounds begin
            U.density_centers[i]      = rhs_density[i]
            U.momentum_centers[i]     = rhs_momentum[i]
            U.total_energy_centers[i] = rhs_energy[i]
        end
    end

end



function LaxFriedrichs_Step_FVS!(W, U::ConservativeVariables, F::FluxVariables,
                                 dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64;
                                 hatF_density = nothing, hatF_momentum = nothing, hatF_total_energy = nothing,
                                 max_wave_speed_fn = nothing)

    # --- prepare / validate scratch arrays (interface indexed 0..total_zones in your layout; we allocate total_zones+1) ---
    if hatF_density === nothing
        hatF_density = zeros(Float64, total_zones+1)
    else
        @assert length(hatF_density) == total_zones+1
    end

    if hatF_momentum === nothing
        hatF_momentum = zeros(Float64, total_zones+1)
    else
        @assert length(hatF_momentum) == total_zones+1
    end

    if hatF_total_energy === nothing
        hatF_total_energy = zeros(Float64, total_zones+1)
    else
        @assert length(hatF_total_energy) == total_zones+1
    end

    # --- build interface fluxes ---
    # We treat `iface` as the interface between cell `iface` and `iface+1`.
    # Valid interface indices: ghost_zones : total_zones - ghost_zones
    Threads.@threads for iface in ghost_zones:(total_zones - ghost_zones)
        @inbounds begin
            # wave speed at this interface
            a_local = if max_wave_speed_fn !== nothing
                max_wave_speed_fn(W, U, iface)
            else
                # fall-back: use magnitudes of center velocities at the two adjacent cells
                max(abs(W.velocity_centers[iface]), abs(W.velocity_centers[iface+1]))
            end

            # density flux at interface `iface`
            hatF_density[iface] =
                0.5*(F.density_flux[iface] + F.density_flux[iface+1]) -
                0.5*a_local*(U.density_centers[iface+1] - U.density_centers[iface])

            # momentum flux
            hatF_momentum[iface] =
                0.5*(F.momentum_flux[iface] + F.momentum_flux[iface+1]) -
                0.5*a_local*(U.momentum_centers[iface+1] - U.momentum_centers[iface])

            # total energy flux
            hatF_total_energy[iface] =
                0.5*(F.total_energy_flux[iface] + F.total_energy_flux[iface+1]) -
                0.5*a_local*(U.total_energy_centers[iface+1] - U.total_energy_centers[iface])
        end
    end

    # --- conservative update on cell centers ---
    # cell centers indices: (ghost_zones+1) : (total_zones - ghost_zones)
    Threads.@threads for i in (ghost_zones+1):(total_zones - ghost_zones)
        @inbounds begin
            U.density_centers[i]      -= (dt/spacing) * (hatF_density[i]      - hatF_density[i-1])
            U.momentum_centers[i]     -= (dt/spacing) * (hatF_momentum[i]     - hatF_momentum[i-1])
            U.total_energy_centers[i] -= (dt/spacing) * (hatF_total_energy[i] - hatF_total_energy[i-1])
        end
    end

end