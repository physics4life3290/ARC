




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

function LaxFriedrichs_Consv_Step!(W::PrimitiveVariables, U::ConservativeVariables, F::FluxVariables, dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64, mode::Symbol, features::Vector{Symbol}, cfl)
    
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
                    # Compute numerical fluxes at interfaces
                    f_density_iphalf = 0.5*(F.density_flux[i] + F.density_flux[i+1]) -
                                    0.5*(spacing/dt)*(U.density_centers[i+1] - U.density_centers[i])
                    f_momentum_iphalf = 0.5*(F.momentum_flux[i] + F.momentum_flux[i+1]) -
                                        0.5*(spacing/dt)*(U.momentum_centers[i+1] - U.momentum_centers[i])
                    f_energy_iphalf   = 0.5*(F.total_energy_flux[i] + F.total_energy_flux[i+1]) -
                                        0.5*(spacing/dt)*(U.total_energy_centers[i+1] - U.total_energy_centers[i])

                    f_density_iminushalf = 0.5*(F.density_flux[i-1] + F.density_flux[i]) -
                                        0.5*(spacing/dt)*(U.density_centers[i] - U.density_centers[i-1])
                    f_momentum_iminushalf = 0.5*(F.momentum_flux[i-1] + F.momentum_flux[i]) -
                                            0.5*(spacing/dt)*(U.momentum_centers[i] - U.momentum_centers[i-1])
                    f_energy_iminushalf   = 0.5*(F.total_energy_flux[i-1] + F.total_energy_flux[i]) -
                                            0.5*(spacing/dt)*(U.total_energy_centers[i] - U.total_energy_centers[i-1])

                    # Update conservative variables
                    U.density_centers[i]      -= (dt/spacing)*(f_density_iphalf   - f_density_iminushalf)
                    U.momentum_centers[i]     -= (dt/spacing)*(f_momentum_iphalf  - f_momentum_iminushalf)
                    U.total_energy_centers[i] -= (dt/spacing)*(f_energy_iphalf    - f_energy_iminushalf)
                end
            end


        elseif mode == :Parallel
        
            Threads.@threads for i in ghost_zones+1:total_zones-ghost_zones
                @inbounds begin
                    # Compute numerical fluxes at interfaces
                    f_density_iphalf = 0.5*(F.density_flux[i] + F.density_flux[i+1]) -
                                    0.5*(spacing/dt)*(U.density_centers[i+1] - U.density_centers[i])
                    f_momentum_iphalf = 0.5*(F.momentum_flux[i] + F.momentum_flux[i+1]) -
                                        0.5*(spacing/dt)*(U.momentum_centers[i+1] - U.momentum_centers[i])
                    f_energy_iphalf   = 0.5*(F.total_energy_flux[i] + F.total_energy_flux[i+1]) -
                                        0.5*(spacing/dt)*(U.total_energy_centers[i+1] - U.total_energy_centers[i])

                    f_density_iminushalf = 0.5*(F.density_flux[i-1] + F.density_flux[i]) -
                                        0.5*(spacing/dt)*(U.density_centers[i] - U.density_centers[i-1])
                    f_momentum_iminushalf = 0.5*(F.momentum_flux[i-1] + F.momentum_flux[i]) -
                                            0.5*(spacing/dt)*(U.momentum_centers[i] - U.momentum_centers[i-1])
                    f_energy_iminushalf   = 0.5*(F.total_energy_flux[i-1] + F.total_energy_flux[i]) -
                                            0.5*(spacing/dt)*(U.total_energy_centers[i] - U.total_energy_centers[i-1])

                    # Update conservative variables
                    U.density_centers[i]      -= (dt/spacing)*(f_density_iphalf   - f_density_iminushalf)
                    U.momentum_centers[i]     -= (dt/spacing)*(f_momentum_iphalf  - f_momentum_iminushalf)
                    U.total_energy_centers[i] -= (dt/spacing)*(f_energy_iphalf    - f_energy_iminushalf)
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

function LaxFriedrichs_Visc_Step!(
    W::PrimitiveVariables,
    U::ConservativeVariables,
    F::FluxVariables,
    dt::Float64,
    ghost_zones::Int,
    total_zones::Int,
    spacing::Float64,
    mode::Symbol,
    features::Vector{Symbol},
    cfl;
    alpha::Float64 = 1.0
)
    LaxFriedrichs_Log = nothing
    U_old = deepcopy(U)

    # Calculate fluxes
    CalculateFlux!(W, U, F)

    if :Debug ∉ features
        if :Verbose ∈ features
            LaxFriedrichs_Log = open("LaxFriedrichs_Log.txt", "a")
            write_solver_input(LaxFriedrichs_Log, W, U_old, F)
        end

        # artificial viscosity coefficient
        ϵ = alpha * spacing / dt

        if mode == :Standard
            for i in ghost_zones+1:total_zones-ghost_zones
                @inbounds begin
                    # Numerical fluxes at interfaces
                    fρ_iphalf = 0.5*(F.density_flux[i] + F.density_flux[i+1]) -
                                0.5*ϵ*(U.density_centers[i+1] - U.density_centers[i])
                    fmu_iphalf = 0.5*(F.momentum_flux[i] + F.momentum_flux[i+1]) -
                                 0.5*ϵ*(U.momentum_centers[i+1] - U.momentum_centers[i])
                    fE_iphalf  = 0.5*(F.total_energy_flux[i] + F.total_energy_flux[i+1]) -
                                 0.5*ϵ*(U.total_energy_centers[i+1] - U.total_energy_centers[i])

                    fρ_iminushalf = 0.5*(F.density_flux[i-1] + F.density_flux[i]) -
                                    0.5*ϵ*(U.density_centers[i] - U.density_centers[i-1])
                    fmu_iminushalf = 0.5*(F.momentum_flux[i-1] + F.momentum_flux[i]) -
                                     0.5*ϵ*(U.momentum_centers[i] - U.momentum_centers[i-1])
                    fE_iminushalf  = 0.5*(F.total_energy_flux[i-1] + F.total_energy_flux[i]) -
                                     0.5*ϵ*(U.total_energy_centers[i] - U.total_energy_centers[i-1])

                    # Update conservative variables
                    U.density_centers[i]      -= (dt/spacing)*(fρ_iphalf  - fρ_iminushalf)
                    U.momentum_centers[i]     -= (dt/spacing)*(fmu_iphalf - fmu_iminushalf)
                    U.total_energy_centers[i] -= (dt/spacing)*(fE_iphalf  - fE_iminushalf)
                end
            end

        elseif mode == :Parallel
            Threads.@threads for i in ghost_zones+1:total_zones-ghost_zones
                @inbounds begin
                    # Numerical fluxes at interfaces
                    fρ_iphalf = 0.5*(F.density_flux[i] + F.density_flux[i+1]) -
                                0.5*ϵ*(U.density_centers[i+1] - U.density_centers[i])
                    fmu_iphalf = 0.5*(F.momentum_flux[i] + F.momentum_flux[i+1]) -
                                 0.5*ϵ*(U.momentum_centers[i+1] - U.momentum_centers[i])
                    fE_iphalf  = 0.5*(F.total_energy_flux[i] + F.total_energy_flux[i+1]) -
                                 0.5*ϵ*(U.total_energy_centers[i+1] - U.total_energy_centers[i])

                    fρ_iminushalf = 0.5*(F.density_flux[i-1] + F.density_flux[i]) -
                                    0.5*ϵ*(U.density_centers[i] - U.density_centers[i-1])
                    fmu_iminushalf = 0.5*(F.momentum_flux[i-1] + F.momentum_flux[i]) -
                                     0.5*ϵ*(U.momentum_centers[i] - U.momentum_centers[i-1])
                    fE_iminushalf  = 0.5*(F.total_energy_flux[i-1] + F.total_energy_flux[i]) -
                                     0.5*ϵ*(U.total_energy_centers[i] - U.total_energy_centers[i-1])

                    # Update conservative variables
                    U.density_centers[i]      -= (dt/spacing)*(fρ_iphalf  - fρ_iminushalf)
                    U.momentum_centers[i]     -= (dt/spacing)*(fmu_iphalf - fmu_iminushalf)
                    U.total_energy_centers[i] -= (dt/spacing)*(fE_iphalf  - fE_iminushalf)
                end
            end

        elseif mode == :GPU
            @warn "GPU mode not yet implemented for Lax–Friedrichs. Defaulting to Parallel..."
            LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)

        elseif mode == :HPC
            @warn "HPC mode not yet implemented for Lax–Friedrichs. Defaulting to Parallel..."
            LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, cfl)
        end

        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, cfl, spacing, dt; logfile=LaxFriedrichs_Log)
            write_solver_output(LaxFriedrichs_Log, W, U, F)
            close(LaxFriedrichs_Log)
        end
    else
        LaxFriedrichs_Step_Debug!(W, U, F, dt, ghost_zones, total_zones, spacing, mode, features, cfl)
    end
end