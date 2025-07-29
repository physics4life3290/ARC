




function Construct1DShockTubePrimitives(_grid::CartesianGrid1D, UserInput)

    wall_positions = sort(UserInput.secondary_input.wall_positions)
    wall_iters = length(wall_positions)
    iters = _grid.xcoord.total_zones
    W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), nothing, nothing, nothing, nothing)   

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
        
        # Now safely index sod_states using region_idx (1-based)
        ρ = UserInput.secondary_input.sod_states[region_idx].ρ
        u = UserInput.secondary_input.sod_states[region_idx].u
        P = UserInput.secondary_input.sod_states[region_idx].P
        γ = UserInput.secondary_input.γ

        W.density_centers[i] = ρ
        W.velocity_centers[i] = u
        W.pressure_centers[i] = P
        W.internal_energy_centers[i] = P / ((γ - 1) * ρ)
    end

    return W
end