

function Construct1DBlastWavePrimitives(_grid::CartesianGrid1D, user_input::UserInput)
    blast_input = user_input.Secondary_Input  # Should be BlastWaveInput
    
    wall_positions = sort(blast_input.wall_positions)
    wall_iters = length(wall_positions)
    iters = _grid.xcoord.total_zones

    # Use the states from parameters.states, which define ambient states
    states = blast_input.parameters.states
    blasts = blast_input.parameters.blasts
    γ = blast_input.gamma

    if user_input.Primary_Input.solver == :FTCS || user_input.Primary_Input.solver == :LaxFriedrichs || user_input.Primary_Input.solver == :Richtmyer
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), nothing, nothing, nothing, nothing)   
    elseif user_input.Primary_Input.solver == :GodunovScheme 
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters)) 
    end

    @inbounds for i in 1:iters
        x = _grid.xcoord.all_centers[i]

        # Determine region index according to walls
        region_idx = 0
        if x <= wall_positions[1]
            region_idx = 1
        elseif x > wall_positions[end]
            region_idx = wall_iters + 1
        else
            for j in 1:wall_iters-1
                if wall_positions[j] < x <= wall_positions[j+1]
                    region_idx = j + 1
                    break
                end
            end
            if region_idx == 0
                error("Unable to assign region index to x = $x")
            end
        end

        ρ = states[region_idx].density_centers
        u = 0.0  # Initial velocity zero
        P = states[region_idx].pressure_centers

        # Add blast energy pressure boost if inside blast radius
        for blast in blasts
            dx = abs(x - blast.center)
            if dx <= blast.radius
                # Distribute blast energy as a simple top-hat pressure increase
                added_e = blast.energy / (2 * blast.radius)  # energy per unit length
                P += (γ - 1) * added_e
            end
        end

        e_int = P / ((γ - 1) * ρ)

        W.density_centers[i] = ρ
        W.velocity_centers[i] = u
        W.pressure_centers[i] = P
        W.internal_energy_centers[i] = e_int
    end

    return W
end
