

function get_BlastWave_states_structured(num_states::Int)
    states = Vector{PrimitiveVariables}(undef, num_states)  

    println("Enter the ambient density (ρ) and pressure (P) for each of the $num_states states.")
    for i in 1:num_states
        while true
            print("State $i → ρ, P: ")
            input = readline()
            values = try
                parse.(Float64, split(input))
            catch
                println("Invalid input. Please enter 2 space-separated numbers.")
                continue
            end

            if length(values) != 2
                println("Please enter exactly 2 values.")
                continue
            end

            states[i] = PrimitiveVariables(values[1], nothing, values[2], nothing, nothing, nothing, nothing, nothing)
            break
        end
    end

    return states
end



function _1DBlastWaveUserInput(grid_input)

    println("Please input how many states are in the box...")
    num_of_states = parse(Int, readline())

    wall_positions = get_wall_positions(num_of_states, grid_input.coord_min, grid_input.coord_max)
    
    states = get_BlastWave_states_structured(num_of_states)

    println("How many blast regions do you want?")
    num_blasts = parse(Int, readline())

    blasts = BlastRegion[]
    for i in 1:num_blasts
        println("Blast $i center position:")
        center = parse(Float64, readline())

        println("Blast $i radius (width of energy deposition):")
        radius = parse(Float64, readline())

        println("Blast $i injected energy:")
        energy = parse(Float64, readline())

        # Ensure blasts are inside the domain
        if center - radius <= grid_input.coord_min || center + radius >= grid_input.coord_max
            println("⚠ Warning: Blast $i is partially or fully outside the domain.")
        end

        push!(blasts, BlastRegion(center, radius, energy))
    end
    
    println("Please input an Adiabatic Constant (γ)...")
    γ = parse(Float64, readline())

    blast_params = BlastParameters(states, blasts)
    
    return BlastWaveInput(blast_params, wall_positions, states, γ)
end
