




function collect_primary_input()
    # Options
    problems = [:ShockTube, :BlastWave, :Custom]
    solvers = [:FTCS, :LaxFriedrichs, :Richtmyer, :GodunovScheme]
    dimensions = ["1D", "2D", "3D"]
    coords = [:Cartesian, :Spherical, :Cylindrical]
    modes = [:Standard, :Parallel, :GPU, :HPC]
    features = [:Animate, :Benchmark, :Debug, :Verbose, :None]
    boundaries = [:Reflecting, :Periodic, :Outflow, :None]
            

    println("Please select the following options to set up your simulation:")
    problem_choice = prompt_choice("Choose from one of the example problems:", problems)

    if problem_choice == :Custom 
        println("When prompted you will need to provide input for density and/or pressure")
        @error "Custom problem setup is not yet implemented. Please choose another problem."
    end

    mode_choice = prompt_choice("Choose a performance mode:", modes)
    feature_choice = prompt_multiple_choices("Choose additional features for solvers (you can select multiple, separated by commas):", features)
    dimension_choice = prompt_choice_string("Choose the dimension of your simulation:", dimensions)
    println("The dimension choice is:", dimension_choice)
    coord_choice = prompt_choice("Choose the coordinate system of your grid:", coords)
    solver_choice = prompt_choice("Choose a numerical solver:", solvers)
    boundary_choice = prompt_choice("Choose boundary conditions:", boundaries)
    println("Enter the name of the output file without the extension...")
    output_filename = readline()

    return PrimaryInput(problem_choice, mode_choice, feature_choice, dimension_choice, coord_choice, solver_choice, boundary_choice, output_filename)
end