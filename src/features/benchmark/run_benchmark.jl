




function run_benchmark(user_input, dens_bench, vel_bench, press_bench)
    if :Benchmark ∈ user_input.Primary_Input.features
        dens_bench .= W.density_centers
        vel_bench .= W.velocity_centers
        press_bench .= W.pressure_centers 

        if user_input.Primary_Input.problem == :ShockTube

            if user_input.Primary_Input.dimension == 1
                # Problem Specific
                pressures = [state.pressure_centers for state in user_input.Secondary_Input.states]
                densities = [state.density_centers for state in user_input.Secondary_Input.states]

                c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
                dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

                W = Construct1DShockTubePrimitives(_grid, user_input)

            end

        elseif user_input.Primary_Input.problem == :BlastWave

            if user_input.Primary_Input.dimension == 1
                pressures = [state.pressure_centers for state in user_input.Secondary_Input.ambient_states]
                densities = [state.density_centers for state in user_input.Secondary_Input.ambient_states]

                c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
                dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

                W = Construct1DBlastWavePrimitives(_grid, user_input)
            end

        end

        U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, 
            W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (user_input.Secondary_Input.gamma - 1)) .+ W.velocity_centers.^2 ./ 2), 
            zeros(length(W.density_centers)), zeros(length(W.density_centers)), zeros(length(W.density_centers)))
        
        F = FluxVariables(zeros(length(U.density_centers)), 
                zeros(length(U.density_centers)),
                zeros(length(U.density_centers)))

        

        while t < t_final

            counter += 1

            ρ, u, p = Riemann_Benchmark(_grid.coord1.centers, _grid.coord1.centers, W.density_centers, W.velocity_centers, W.pressure_centers, user_input.Secondary_Input.gamma, t)
            
            F.density_flux = ρ .* u
            F.momentum_flux = ρ .* u .^ 2 .+ p
            F.total_energy_flux = u .* (p ./ (user_input.Secondary_Input.gamma - 1) .+ 0.5 .* ρ .* u .^ 2 .+ p)

            source_terms = nothing
            if user_input.Primary_Input.coordinate_system == :cylindrical
                source_terms = (1/_grid.coord1.all_centers) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
            elseif user_input.Primary_Input.coordinate_system == :spherical
                source_terms = (2/(_grid.coord1.all_centers)) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
            end

            W.density_centers .= ρ
            W.velocity_centers .= u
            W.pressure_centers .= p
            W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2

            # This Clamping is for oscillatory schemes like FTCS and Richtmyer
            Threads.@threads for i in 1:_grid.coord1.total_zones
                @inbounds begin
                    if W.pressure_centers[i] < 0.0
                        W.pressure_centers[i] = 0.0
                    elseif W.density_centers[i] < 0.0
                        W.density_centers[i] = 0.0
                    elseif W.internal_energy_centers[i] < 0.0
                        W.internal_energy_centers[i] = 0.0   
                    end
                end
            end

            

            c = sqrt.(user_input.Secondary_Input.gamma .* W.pressure_centers ./ W.density_centers)
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / maximum(abs.(W.velocity_centers) .+ c)
            dt = min(dt, t_final - t)
            t += dt
        end
        dens_exact, vel_exact, press_exact = W.density_centers, W.velocity_centers, W.pressure_centers

        dens_norms = compute_norms(dens_bench, dens_exact)
        vel_norms = compute_norms(vel_bench, vel_exact)
        press_norms = compute_norms(press_bench, press_exact)

        println("Comparing the results from the $(user_input.Primary_Input.solver) solver to the exact Riemann solution: ")
        println("Density norms: ", dens_norms)
        println("Velocity norms: ", vel_norms)
        println("Pressure norms: ", press_norms)

    end

end