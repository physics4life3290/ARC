




PROBLEMS   = [:ShockTube, :BlastWave, :Custom]
SOLVERS    = [:FTCS, :LaxFriedrichs, :Richtmyer, :GodunovScheme]
DIMENSIONS = ["1D", "2D", "3D"]
COORDS     = [:Cartesian, :Spherical, :Cylindrical]
MODES      = [:Standard, :Parallel, :GPU, :HPC]
FEATURES   = [:Animate, :Benchmark, :Debug, :Verbose, :None]
BOUNDARIES = [:Reflecting, :Periodic, :Outflow, :None]


function collect_primary_input()
    println("Please select the following options to set up your simulation:")

    problem_choice   = prompt_choice("Choose from one of the example problems:", PROBLEMS)
    if problem_choice == :Custom
        println("When prompted you will need to provide input for density and/or pressure")
        @error "Custom problem setup is not yet implemented. Please choose another problem."
    end

    mode_choice      = prompt_choice("Choose a performance mode:", MODES)
    feature_choice   = prompt_multiple_choices("Choose additional features (comma separated):", FEATURES)
    dimension_choice = prompt_choice_string("Choose the dimension of your simulation:", DIMENSIONS)
    coord_choice     = prompt_choice("Choose the coordinate system of your grid:", COORDS)
    solver_choice    = prompt_choice("Choose a numerical solver:", SOLVERS)

    boundary_choice  = (:Benchmark in feature_choice) ? :None :
                      prompt_choice("Choose boundary conditions:", BOUNDARIES)

    println("Enter the name of the output file without the extension...")
    output_filename = readline()

    return PrimaryInput(problem_choice, mode_choice, feature_choice,
                        dimension_choice, coord_choice, solver_choice,
                        boundary_choice, output_filename)
end
