




function FTCS_Step_Debug!(W, U::ConservativeVariables,
                    F::FluxVariables,
                    dt::Float64,
                    ghost_zones::Int,
                    total_zones::Int,
                    spacing::Float64, 
                    mode::Symbol, 
                    features::Vector{Symbol}, 
                    CFL::Float64
                    )

    println("You are debugging the FTCS step function...")
    # 1) copy old solution
    try
        U_old = deepcopy(U)
        CalculateFlux!(W, U_old, F)

        if :Verbose ∈ features
            println("You are using the FTCS method. This method is Unconditionally unstable for the hyperbolic conservation laws.\n It's sole purpose is educational.")
            println("Opening log file for FTCS...")
            FTCS_Log = open("FTCS_Log.txt", "a")
            write_solver_input(FTCS_Log, W, U_old, F)
        end

        # 2) compute new solution from old    
        if mode == :Standard

            for i in (ghost_zones+1):(total_zones-ghost_zones)

                @inbounds begin
                    U.density_centers[i]      = U_old.density_centers[i]      - (dt/(2*spacing))*(F.density_flux[i+1]  - F.density_flux[i-1])
                    U.momentum_centers[i]     = U_old.momentum_centers[i]     - (dt/(2*spacing))*(F.momentum_flux[i+1]  - F.momentum_flux[i-1])
                    U.total_energy_centers[i] = U_old.total_energy_centers[i] - (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
                end

            end

        elseif mode == :Parallel

            Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
                @inbounds begin
                    U.density_centers[i]      = U_old.density_centers[i]      - (dt/(2*spacing))*(F.density_flux[i+1]  - F.density_flux[i-1])
                    U.momentum_centers[i]     = U_old.momentum_centers[i]     - (dt/(2*spacing))*(F.momentum_flux[i+1]  - F.momentum_flux[i-1])
                    U.total_energy_centers[i] = U_old.total_energy_centers[i] - (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
                end
            end

        elseif mode == :GPU

            @warn "GPU mode not yet implemented for FTCS... Release coming in v0.2... Defaulting to Parallel..."
            FTCS_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, CFL)

        elseif mode == :HPC

            @warn "HPC mode not yet implemented for FTCS... Release coming in v0.5...Defaulting to Parallel..."
            FTCS_Step!(W, U, F, dt, ghost_zones, total_zones, spacing, :Parallel, features, CFL)

        end

        if :Verbose ∈ features
            solver_diagnostics(U, U_old, F, CFL, spacing, dt; logfile=FTCS_Log)
            write_solver_output(FTCS_Log, W, U, F)
            close(FTCS_Log)
        end
        
    catch e 
        println("The incoming data is 
        Density: $(U_old.density_centers)
        Length of Density: $(length(U_old.density_centers))
        Momentum: $(U_old.momentum_centers)
        Length of Momentum: $(length(U_old.momentum_centers))
        Total Energy: $(U_old.total_energy_centers)
        Length of Total Energy: $(length(U_old.total_energy_centers))")

        @error "An error occured during FTCS step: $e"
    end
end