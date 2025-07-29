

# Prompt user
function run_introduction()
    println(
    "    #====================#
    #  The CΩδΞχ Trials  #
    #====================#")
    println("This is the test suite for ARC...")
    println("Press enter to continue...")
    readline()
end



function initiate_UI()
   
    # Options
    #problems = ["Sod Shock/Riemann Type Problem", "Sedov Blast Wave", "Shu-Osher Problem"]
    problems = [:ShockTube, :BlastWave, :ShuOsher]
    #solvers = ["FTCS", "Lax-Friedrichs", "Richtmyer", "Godunov's Method", "MUSCL", "PPM"]
    solvers = [:FTCS, :LaxFriedrichs, :Richtmyer, :GodunovsScheme, :MUSCL, :PPM]
    dimensions = ["1D", "2D"]
    #coords = ["Cartesian", "Spherical", "Cylindrical"]
    coords = [:Cartesian, :Spherical, :Cylindrical]
    #modes = ["Standard (HDF5 output)", "Verbose (HDF5 and general .txt output)", "Debug (In depth .txt output)", "Animate (GIF and HDF5 output)", "Benchmark ()"]
    modes = [:Standard, :Verbose, :Debug, :Animate, :Benchmark]
    #boundaries = ["Reflecting", "Periodic", "Outflow", "No Boundary Condition"]
    boundaries = [:Reflecting, :Periodic, :Outflow, :None]

    run_introduction()

    println("Please select the following options to set up your simulation:")
    problem_choice = prompt_choice("Choose a problem:", problems)
    mode_choice = prompt_choice("Choose a mode:", modes)
    dimension_choice = prompt_choice_string("Choose the dimension:", dimensions)
    coord_choice = prompt_choice("Choose the coordinate system:", coords)
    solver_choice = prompt_choice("Choose a solver:", solvers)
    boundary_choice = prompt_choice("Choose boundary conditions:", boundaries)
    println("Enter the name of the output file without the extension...")
    output_filename = readline()

    primary_input = (problem = problem_choice, mode = mode_choice, filename = output_filename, dimension = dimension_choice, coord_sys = coord_choice, solver = solver_choice, boundary_condition = boundary_choice)
    
    
    if primary_input.problem == :ShockTube
        if primary_input.dimension == 1
            if primary_input.coord_sys == :Cartesian
                println("You have selected 1D Sod Shock Problem in Cartesian Coordinates...")
                secondary_input = (domain = nothing, grid_center = nothing, zones = nothing, 
                                    ghost_zones = nothing, x_min = nothing, x_max = nothing, 
                                    wall_positions = nothing, sod_states = nothing, γ = nothing, cfl = nothing, t_final=nothing)
                secondary_input = _1DSODUserInput(primary_input)
                UserInput = (primary_input=primary_input, secondary_input=secondary_input)
            end
        end
    end
    
    return UserInput
end