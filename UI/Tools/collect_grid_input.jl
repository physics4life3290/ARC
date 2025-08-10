





function collect_grid_input(primary_input)
        if primary_input.dimension == 1
            if primary_input.coordinate_system == :Cartesian
                println("Please choose the parameters for the 1D Cartesian Grid...")
                println("Please input the domain length...")
                domain_length = parse(Float64, readline())
                println("Please input the point the grid is centered about...")
                grid_center = parse(Float64, readline())
                println("Please input the number of zones in the grid...")
                zones = parse(Int, readline())
                println("Please input the number of ghost zones...")
                ghost_zones = parse(Int, readline())
                total_zones = zones + 2 * ghost_zones
                println("The grid will contain $(total_zones) total zones...")
                x_min = -domain_length/2 + grid_center
                x_max = domain_length/2 + grid_center
                return GridInput(domain_length, grid_center, zones, ghost_zones, total_zones, x_min, x_max)
            end
        end
    end