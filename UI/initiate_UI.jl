




function initiate_UI()
   
    # Options
    problems = ["Sod Shock Tube", "Sedov Blast Wave", "Shu-Osher Problem", "Double Blast Wave Interaction", "Double Mach Reflection"]
    solvers = ["FTCS", "Lax-Friedrichs", "Richtmyer"]
    dimensions = ["1D", "2D"]
    coords = ["Cartesian", "Spherical"]
    boundaries = ["Reflecting", "Periodic", "Outflow", "No Boundary Condition"]

    # Prompt user
    println("Welcome to The Codex Trials...")
    println("This is the test suite for ARC...")
    problem_choice = prompt_choice("Choose a problem:", problems)
    dimension_choice = prompt_choice("Choose the dimension:", dimensions)
    coord_choice = prompt_choice("Choose the coordinate system:", coords)
    num_zones = prompt_number("Enter the number of zones in the grid:", Int)
    num_ghost_zones = prompt_number("Enter the number of ghost zones:", Int)
    boundary_choice = prompt_choice("Choose a boundary condition:", boundaries)
    solver_choice = prompt_choice("Choose a solver:", solvers)
    if solver_choice == 1
        println("You have chosen the FTCS solver.")
        println("Note: FTCS is conditionally stable and will not converge...")
        println("Its purpose is educational, not practical.")
    end
    cfl = prompt_number("Enter the CFL number:", Float64)
    if cfl <= 0
        println("CFL number must be positive. Defaulting to 0.5.")
        cfl = 0.5
    elseif cfl > 1
        println("CFL number should typically be less than or equal to 1 as it may lead to instability in some cases. Proceeding with the provided value.")
    end
    t_final = prompt_number("Enter the final time for the simulation:", Float64)
    print("Enter the output filename (without extension): ")
    output_filename = readline()

    # Now that user input is collected for the general problem, we can handle specific problem setup
    if problem_choice == 1
        println("You have chosen the Sod Shock Tube problem.")
        # We want to build a function to set up the Sod Shock Tube problem
    elseif problem_choice == 2
        println("You have chosen the Sedov Blast Wave problem.")
    elseif problem_choice == 3
        println("You have chosen the Shu-Osher Problem.")
    elseif problem_choice == 4
        println("You have chosen the Double Blast Wave Interaction problem.")
    elseif problem_choice == 5
        println("You have chosen the Double Mach Reflection problem.")
        if dimension_choice < 2
            println("Double Mach Reflection is only available in 2D. Defaulting to 2D.")
            dimension_choice = 2
        end
    end

    # Collect parameters
    params = (
        problem = problems[problem_choice],
        solver = solvers[solver_choice],
        dimension = dimensions[dimension_choice],
        coordinate_system = coords[coord_choice],
        boundary_condition = boundaries[boundary_choice],
        num_zones = num_zones,
        num_ghost_zones = num_ghost_zones,
        cfl = cfl,
        t_final = t_final,
        output_filename = output_filename,
    )

    # Print selections
    println("\nYou selected:")
    for (k,v) in pairs(params)
        println("  $k ---> $v")
    end

    return params
end