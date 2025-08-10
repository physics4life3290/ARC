




function collect_primary_input()
    # Options
    problems = [:ShockTube, :BlastWave]
    solvers = [:FTCS, :LaxFriedrichs, :Richtmyer, :GodunovScheme]
    dimensions = ["1D", "2D"]
    coords = [:Cartesian, :Spherical, :Cylindrical]
    modes = [:Standard, :Verbose, :Debug, :Animate, :Benchmark]
    boundaries = [:Reflecting, :Periodic, :Outflow, :None]
            

    println("Please select the following options to set up your simulation:")
    problem_choice = prompt_choice("Choose from one of the example problems:", problems)
    mode_choice = prompt_choice("Choose a performance mode:", modes)
    dimension_choice = prompt_choice_string("Choose the dimension of your simulation:", dimensions)
    coord_choice = prompt_choice("Choose the coordinate system of your grid:", coords)
    solver_choice = prompt_choice("Choose a numerical solver:", solvers)
    boundary_choice = prompt_choice("Choose boundary conditions:", boundaries)
    println("Enter the name of the output file without the extension...")
    output_filename = readline()

    return PrimaryInput(problem_choice, mode_choice, dimension_choice, coord_choice, solver_choice, boundary_choice, output_filename)
end