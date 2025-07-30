
using BenchmarkTools



function run_Shock_Tube_Benchmark(UserInput)
    
    if UserInput.primary_input.dimension == 1
        zones = parse(Int, UserInput.secondary_input.zones)
        ghost_zones = parse(Int, UserInput.secondary_input.ghost_zones)
        grid_center = UserInput.secondary_input.grid_center
        
        init_grid = time()
        _grid = Construct1DCartesian(UserInput.secondary_input.domain, zones, ghost_zones, grid_center, "cm")
        final_grid = time()

        init_time_step = time()
        pressures = [state.P for state in UserInput.secondary_input.sod_states]
        densities = [state.ρ for state in UserInput.secondary_input.sod_states]
        c = sqrt(UserInput.secondary_input.γ * maximum(pressures) / maximum(densities))
        dt = UserInput.secondary_input.cfl * _grid.xcoord.spacing / c
        final_time_step = time()

        init_prim_vars = time()
        W = Construct1DShockTubePrimitives(_grid, UserInput)
        final_prim_vars = time()
        
        init_conserv_vars = time()
        U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (UserInput.secondary_input.γ - 1)) .+ W.velocity_centers.^2 ./ 2), nothing, nothing, nothing)
        final_conserv_vars = time()

        F = FluxVariables(zeros(length(U.density_centers)), 
                zeros(length(U.density_centers)),
                zeros(length(U.density_centers)))

        counter = 0
        t = 0.0
        t_final = parse(Float64, UserInput.secondary_input.t_final)

        init_evolution = time()
        while t < t_final
            t += dt
            counter += 1
            apply_boundary_conditions(UserInput, U, _grid)

            F.density_flux .= U.momentum_centers 
            F.momentum_flux .= U.momentum_centers .* W.velocity_centers .+ W.pressure_centers
            F.total_energy_flux .= W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)

            if UserInput.primary_input.solver == :FTCS
                FTCS_Step!(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
            elseif UserInput.primary_input.solver == :LaxFriedrichs
                LaxFriedrichs_Step(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
            elseif UserInput.primary_input.solver == :Richtmyer
                RichtmyerStep!(U, F, _grid, UserInput,dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing, UserInput.secondary_input.γ)
            else 
                println("Defaulting to Lax until Scheme requested is supported...")
                LaxFriedrichs_Step(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
            end

            W.density_centers .= U.density_centers
            W.velocity_centers .= U.momentum_centers ./ U.density_centers
            W.pressure_centers .= (UserInput.secondary_input.γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
            W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
            

            # This Clamping is for oscillatory schemes like FTCS and Richtmyer
            @inbounds for i in 1:_grid.xcoord.total_zones
                if W.pressure_centers[i] < 0.0
                    W.pressure_centers[i] = 0.0
                elseif W.density_centers[i] < 0.0
                    W.density_centers[i] = 0.0
                elseif W.internal_energy_centers[i] < 0.0
                    W.internal_energy_centers[i] = 0.0   
                end
            end

            cs = sqrt.(UserInput.secondary_input.γ .* W.pressure_centers ./ W.density_centers)
            wavespeed = abs.(W.velocity_centers) .+ cs
            maxspeed = maximum(wavespeed)
            dt = UserInput.secondary_input.cfl * _grid.xcoord.spacing / max(maxspeed, 1e-6)
            
        end
        final_evolution = time()

        benchmark_log = open("Benchmark_Results.txt", "w")
        write(benchmark_log, "The time to create the grid is: $(final_grid - init_grid)
The time to calculate the time step is: $(final_time_step - init_time_step)
The time to configure the Primitive and Conservative vars are: $(final_prim_vars - init_prim_vars), $(final_conserv_vars - init_conserv_vars)
The time it takes for the whole evolution is: $(final_evolution - init_evolution)")
            close(benchmark_log)
    end
end