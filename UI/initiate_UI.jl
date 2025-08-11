

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

   run_introduction()

    primary_input = collect_primary_input()
    
    grid_input = collect_grid_input(primary_input)

    solver_input = collect_solver_input(primary_input)

    secondary_input = collect_secondary_input(primary_input, grid_input)
    
    user_input = UserInput(primary_input, grid_input, solver_input, secondary_input)
    return user_input
end