




function Construct1DShockTubePrimitives(_grid::CartesianGrid1D, user_input)

    wall_positions = sort(user_input.Secondary_Input.wall_positions)
    wall_iters = length(wall_positions)
    iters = _grid.xcoord.total_zones
    if user_input.Primary_Input.solver == :FTCS || user_input.Primary_Input.solver == :LaxFriedrichs || user_input.Primary_Input.solver == :Richtmyer
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), nothing, nothing, nothing, nothing)   
    elseif user_input.Primary_Input.solver == :GodunovScheme 
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters)) 
    end
        @inbounds for i in 1:iters
        x = _grid.xcoord.all_centers[i]

        region_idx = 0

        if x <= wall_positions[1]
            region_idx = 1
        elseif x > wall_positions[end]
            region_idx = wall_iters + 1
        else
            @inbounds for j in 1:wall_iters-1
                if wall_positions[j] < x <= wall_positions[j+1]
                    region_idx = j + 1
                    break
                end
            end
            # x is between wall_positions[1] and wall_positions[end] but not between any pair?
            if region_idx == 0
                @inbounds for j in 1:wall_iters-1
                    println("DEBUG: x = $x, wall range = ($(wall_positions[j]), $(wall_positions[j+1]))")
                end
                error("Unable to assign region index to x = $x")
            end
        end
        
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

    return W
end