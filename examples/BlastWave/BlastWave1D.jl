

function Construct1DBlastWavePrimitives(_grid::CartesianGrid1D, user_input::UserInput)
    blast_input = user_input.Secondary_Input  # Should be BlastWaveInput
    
    wall_positions = sort(blast_input.wall_positions)
    wall_iters = length(wall_positions)
    iters = _grid.coord1.total_zones

    # Use the states from parameters.states, which define ambient states
    states = blast_input.parameters.states
    blasts = blast_input.parameters.blasts
    γ = blast_input.gamma

    if user_input.Primary_Input.solver == :FTCS || user_input.Primary_Input.solver == :LaxFriedrichs || user_input.Primary_Input.solver == :Richtmyer
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), nothing, nothing, nothing, nothing)   
    elseif user_input.Primary_Input.solver == :GodunovScheme 
        W = PrimitiveVariables(zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters), zeros(iters)) 
    end

    for i in 1:iters
        @inbounds begin 
            x = _grid.coord1.all_centers[i]

            # Determine region index according to walls
            
            region_idx  = get_regions(x, wall_positions, wall_iters)
            ρ = states[region_idx].density_centers
            u = 0.0  # Initial velocity zero
            P = states[region_idx].pressure_centers

            # Add blast energy pressure boost if inside blast radius
            P = get_blasts(blasts)
            e_int = P / ((γ - 1) * ρ)

            W.density_centers[i] = ρ
            W.velocity_centers[i] = u
            W.pressure_centers[i] = P
            W.internal_energy_centers[i] = e_int
        end
    end

    return W
end
