




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

struct FluidState
    ρ::Float64
    u::Float64
    P::Float64
end

function get_states_structured(num_states::Int)
    states = Vector{FluidState}(undef, num_states)

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

            states[i] = FluidState(values[1], values[2], values[3])
            break
        end
    end

    return states
end

function _1DSODUserInput(primary_input)
    println("Please input the length of the box...")
    domain_length = readline()
    domain_length = parse(Float64, domain_length)
    println("Please input grid center...")
    grid_center = readline()
    grid_center = parse(Float64, grid_center)
    println("Please input the number of grid points...")
    zones = readline()
    println("Please input the number of ghost zones...")
    ghost_zones = readline()
    println("Please input how many states are in the box...")
    num_of_states = readline()
    num_of_states = parse(Int, num_of_states)
    x_min = -1*domain_length/2 + grid_center
    x_max = domain_length/2 + grid_center
    wall_positions = get_wall_positions(num_of_states, x_min, x_max)
    states = get_states_structured(num_of_states)
    println("Please input an Adiabatic Constant...")
    adiabat_const = readline()
    adiabat_const = parse(Float64, adiabat_const)
    println("Please input the Courant Condition...")
    CFL = readline()
    CFL = parse(Float64, CFL)
    println("Please select the final time...")
    final_t = readline()      
    return (domain = domain_length, grid_center = grid_center, zones = zones, 
    ghost_zones = ghost_zones, x_min = x_min, x_max = x_max, wall_positions = wall_positions, 
    sod_states = states, γ = adiabat_const, cfl = CFL, t_final = final_t)
end