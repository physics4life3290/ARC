



function parse_point(str::String)::Tuple{Float64,Float64}
    # Remove parentheses and extra spaces
    s = replace(str, r"[()]" => "")
    s = strip(s)

    # Split on either comma or whitespace
    parts = split(s, r"[\s,]+")

    if length(parts) != 2
        error("Invalid point format. Expected two numbers, got: $str")
    end

    return (parse(Float64, parts[1]), parse(Float64, parts[2]))
end


function collect_grid_input(primary_input)
    if primary_input.dimension == 1
        if primary_input.coordinate_system == :Cartesian
            
            println("Please choose the parameters for the 1D Cartesian Grid...")
            domain_length = prompt_domain_length()
            
            println("Please input the point the grid is centered about...")
            grid_center = parse(Float64, readline())
            
            zones = prompt_zones()
            
            if primary_input.solver == :FTCS || primary_input.solver == :LaxFriedrichs || primary_input.solver == :Richtmyer
                println("Please input the number of ghost zones (must be at least one)...")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            elseif primary_input.solver == :GodunovScheme 
                println("Please input the number of ghost zones (must be at least one for First Order, two for MUSCL, and three for PPM)")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            end

            println("The grid will contain $(total_zones) total zones...")
            x_min = -domain_length/2 + grid_center
            x_max = domain_length/2 + grid_center
            
            return GridInput((domain_length,), (grid_center,), (zones,), ghost_zones, (total_zones,), (x_min,), (x_max,))

        elseif primary_input.coordinate_system == :Spherical || primary_input.coordinate_system == :Cylindrical
            println("Please choose the parameters for the 1D Spherical (or Cylindrical) Grid...")
            println("Please input the domain length...")
            domain_length = parse(Float64, readline())
            grid_center = prompt_grid_min()

            zones = prompt_zones()

            if primary_input.solver == :FTCS || primary_input.solver == :LaxFriedrichs || primary_input.solver == :Richtmyer
                println("Please input the number of ghost zones (must be at least one)...")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            elseif primary_input.solver == :GodunovScheme 
                println("Please input the number of ghost zones (must be at least one for First Order, two for MUSCL, and three for PPM)")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            end

            println("The grid will contain $(total_zones) total zones...")
            r_min = 1E-12 + grid_center
            r_max = domain_length + grid_center
            return GridInput((domain_length, ), (grid_center, ), (zones, ), ghost_zones, (total_zones, ), (r_min, ), (r_max, ))
        end

    elseif primary_input.dimension == 2
        if primary_input.coordinate_system == :Cartesian 
            println("Please choose the parameters for the 2D Cartesian Grid...(Press enter to continue!)")
            readline()
            println("Select the domain for the x coordinate...")
            x_domain_length = prompt_domain_length()
            println("Select the domain for the y coordinate...")
            y_domain_length = prompt_domain_length()

            println("Please input the point the grid is centered about...")
            grid_center = parse_point(readline())

            println("Parsed center: ", grid_center)

            println("Please input the number of x zones...")
            x_zones = prompt_zones()
            println("Please input the number of y zones...")
            y_zones = prompt_zones()

            if primary_input.solver == :FTCS || primary_input.solver == :LaxFriedrichs || primary_input.solver == :Richtmyer
                println("Please input the number of ghost zones (must be at least one)...")
                ghost_zones = parse(Int, readline())
                x_total_zones = (x_zones + 2 * ghost_zones)
                y_total_zones = (y_zones + 2 * ghost_zones)
            elseif primary_input.solver == :GodunovScheme 
                println("Please input the number of ghost zones (must be at least one for First Order, two for MUSCL, and three for PPM)")
                ghost_zones = parse(Int, readline())
                x_total_zones = (x_zones + 2 * ghost_zones)
                y_total_zones = (y_zones + 2 * ghost_zones)
            end

            println("The grid will contain $(x_total_zones*y_total_zones) total zones...")
            x_min = -x_domain_length/2 + grid_center[1]
            x_max = x_domain_length/2 + grid_center[1]
            y_min = -y_domain_length/2 + grid_center[2]
            y_max = y_domain_length/2 + grid_center[2]
                
            return GridInput((x_domain_length, y_domain_length), grid_center, (x_zones, y_zones), ghost_zones, (x_total_zones, y_total_zones), (x_min, y_min), (x_max, y_max))

        elseif primary_input.coordinate_system == :Spherical || primary_input.coordinate_system == :Cylindrical
            println("Please choose the parameters for the 1D Spherical (or Cylindrical) Grid...")
            println("Please input the domain length...")
            domain_length = parse(Float64, readline())
            grid_center = prompt_grid_min()

            zones = prompt_zones()

            if primary_input.solver == :FTCS || primary_input.solver == :LaxFriedrichs || primary_input.solver == :Richtmyer
                println("Please input the number of ghost zones (must be at least one)...")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            elseif primary_input.solver == :GodunovScheme 
                println("Please input the number of ghost zones (must be at least one for First Order, two for MUSCL, and three for PPM)")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
            end

            println("The grid will contain $(total_zones) total zones...")
            r_min = 1E-12 + grid_center
            r_max = domain_length + grid_center
            return GridInput(domain_length, grid_center, zones, ghost_zones, total_zones, r_min, r_max)
        end

    end
end