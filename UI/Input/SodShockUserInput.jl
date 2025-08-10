




function get_wall_positions(num_states::Int, x_min::Float64, x_max::Float64)
    num_walls = num_states - 1
    wall_positions = Float64[]
    
    println("Your domain spans from $x_min to $x_max.")
    println("Enter $num_walls wall positions strictly within this domain (e.g. -0.5 0.0 0.5):")
    
    while true
        input = readline()
        try
            positions = sort(parse.(Float64, split(input)))
            if length(positions) != num_walls
                println("You must enter exactly $num_walls values.")
                continue
            end
            if any(pos -> pos ≤ x_min || pos ≥ x_max, positions)
                println("All wall positions must lie strictly within the domain ($x_min, $x_max).")
                continue
            end
            return positions
        catch e
            println("Invalid input. Please enter $num_walls space-separated numbers.")
        end
    end
end

function get_ShockTube_states_structured(num_states::Int)
    states = Vector{PrimitiveVariables}(undef, num_states)  

    println("Enter density (ρ), velocity (u), and pressure (P) for each of the $num_states states.")
    for i in 1:num_states
        while true
            print("State $i → ρ, u, P: ")
            input = readline()
            values = try
                parse.(Float64, split(input))
            catch
                println("Invalid input. Please enter 3 space-separated numbers.")
                continue
            end

            if length(values) != 3
                println("Please enter exactly 3 values.")
                continue
            end

            states[i] = PrimitiveVariables(values[1], values[2], values[3], nothing, nothing, nothing, nothing, nothing)
            break
        end
    end

    return states
end

function _1DShockTubeUserInput(grid_input)

    println("Please input how many states are in the box...")
    num_of_states = parse(Int, readline())

    wall_positions = get_wall_positions(num_of_states, grid_input.coord_min, grid_input.coord_max)

    states = get_ShockTube_states_structured(num_of_states)

    println("Please input an Adiabatic Constant...")
    adiabat_const = parse(Float64, readline())

    return ShockTubeInput(wall_positions, states, adiabat_const)

end