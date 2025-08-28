





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
                
                return GridInput(domain_length, grid_center, zones, ghost_zones, total_zones, x_min, x_max)

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