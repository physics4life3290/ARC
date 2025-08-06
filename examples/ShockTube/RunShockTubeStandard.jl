




function run_Shock_Tube_Standard(UserInput)
    println("The run_Shock_Tube function has been created and called successfully!")
    h5_filename = UserInput.primary_input.filename
    if UserInput.primary_input.dimension == 1
        zones = parse(Int, UserInput.secondary_input.zones)
        ghost_zones = parse(Int, UserInput.secondary_input.ghost_zones)
        grid_center = UserInput.secondary_input.grid_center
        
        _grid = Construct1DCartesian(UserInput.secondary_input.domain, zones, ghost_zones, grid_center, "cm")
        println("The grid has been successfully created!")
        pressures = [state.P for state in UserInput.secondary_input.sod_states]
        densities = [state.ρ for state in UserInput.secondary_input.sod_states]

        c = sqrt(UserInput.secondary_input.γ * maximum(pressures) / maximum(densities))
        dt = UserInput.secondary_input.cfl * _grid.xcoord.spacing / c

        W = Construct1DShockTubePrimitives(_grid, UserInput)
        
        U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (UserInput.secondary_input.γ - 1)) .+ W.velocity_centers.^2 ./ 2), nothing, nothing, nothing)

        F = FluxVariables(zeros(length(U.density_centers)), 
                zeros(length(U.density_centers)),
                zeros(length(U.density_centers)))

        counter = 0
        t = 0.0
        t_final = parse(Float64, UserInput.secondary_input.t_final)

        while t < t_final
            t += dt
            counter += 1

            apply_boundary_conditions(UserInput, U, _grid)
            
            if UserInput.primary_input.solver == :FTCS
                FTCS_Step!(W, U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
            elseif UserInput.primary_input.solver == :LaxFriedrichs
                LaxFriedrichs_Step(W, U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
            elseif UserInput.primary_input.solver == :Richtmyer
                RichtmyerStep!(W, U, F, _grid, UserInput,dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing, UserInput.secondary_input.γ)
            elseif UserInput.primary_input.solver == :GodunovsScheme
                Godunov_Step!(UserInput, _grid, W, U)
            else 
                println("Defaulting to Lax until Scheme requested is supported...")
                LaxFriedrichs_Step(W, U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
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
           
            if counter % 1 == 0
                groupname = "step_$(counter)"
                println("Saving snapshot to $h5_filename in group $groupname")
                mode = isfile(h5_filename) ? "r+" : "w"
                h5open(h5_filename, mode) do file  # "a" = append mode
                    grp = create_group(file, groupname)

                    grp["x"] = _grid.xcoord.all_centers

                    grp["W/density"]      = W.density_centers
                    grp["W/velocity"]     = W.velocity_centers
                    grp["W/pressure"]     = W.pressure_centers
                    grp["W/int_energy"]   = W.internal_energy_centers

                    grp["U/density"]      = U.density_centers
                    grp["U/momentum"]     = U.momentum_centers
                    grp["U/total_energy"] = U.total_energy_centers

                    grp["F/dens_flux"]    = F.density_flux
                    grp["F/mome_flux"]    = F.momentum_flux
                    grp["F/tot_ener_flux"]= F.total_energy_flux

                    grp["time"] = t
                end
                
            end
            
        end
    end
end