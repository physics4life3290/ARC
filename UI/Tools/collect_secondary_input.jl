




function collect_secondary_input(primary_input, grid_input)
    if primary_input.dimension == 1
        if primary_input.coordinate_system == :Cartesian
            if primary_input.problem == :ShockTube
                println("You have selected 1D Shock Tube Problem in Cartesian Coordinates...")
                secondary_input = (wall_positions = nothing, states = nothing, γ = nothing)
                return _1DShockTubeUserInput(grid_input)
            elseif primary_input.problem == :BlastWave
                println("You have selected 1D Blast Wave Problem in Cartesian Coordinates...")
                secondary_input = (blast_params = nothing, wall_positions = nothing, states = nothing, γ = nothing)
                return _1DBlastWaveUserInput(grid_input)
            end
        end
    end
end