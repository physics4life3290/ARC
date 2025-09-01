




function Construct1DShockTubePrimitives(_grid::Union{CartesianGrid1D, SphericalGrid1D}, user_input)

    wall_positions = sort(user_input.Secondary_Input.wall_positions)
    wall_iters = length(wall_positions)
    iters = _grid.coord1.total_zones
    if user_input.Primary_Input.solver == :FTCS || user_input.Primary_Input.solver == :LaxFriedrichs || user_input.Primary_Input.solver == :Richtmyer
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), nothing, nothing, nothing, nothing)   
    elseif user_input.Primary_Input.solver == :GodunovScheme 
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters)) 
    end
    
    for i in 1:iters
        @inbounds begin 
            x = _grid.coord1.all_centers[i]
            region_idx = get_regions(x, wall_positions, wall_iters)
            # Now safely index states using region_idx (1-based)
            ρ = user_input.Secondary_Input.states[region_idx].density_centers
            u = user_input.Secondary_Input.states[region_idx].velocity_centers
            P = user_input.Secondary_Input.states[region_idx].pressure_centers
            γ = user_input.Secondary_Input.gamma

            W.density_centers[i] = ρ
            W.velocity_centers[i] = u
            W.pressure_centers[i] = P
            W.internal_energy_centers[i] = P / ((γ - 1) * ρ)
        end
    end

    return W
end