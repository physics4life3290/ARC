




function validate_mode_choice(mode_choice, output_filename)
    if mode_choice == 1
        println("You have chosen Fast mode.")
        output_filename *= ".h5" #produces an HDF5 file
    elseif mode_choice == 2
        println("You have chosen Show Details mode.")
        output_filename *= ".txt" #produces a text file with details
    elseif mode_choice == 3
        println("You have chosen Save and Log mode.")
        output_filename_log = output_filename
        output_filename *= ".h5" #produces an HDF5 file
        output_filename_log *= ".txt" #produces a text file with details
    elseif mode_choice == 4
        println("You have chosen Animate mode.")
        output_filename *= ".gif" #produces an animated GIF
    end
    return output_filename
end

function validate_solver_choice(solver_choice)
    if solver_choice == 1
        println("You have chosen the FTCS solver.")
        @warn "Note: FTCS is unconditionally unstable and will not converge..."
        println("Its purpose is educational, not practical.")
    elseif solver_choice == 2
        println("You have chosen the Lax-Friedrichs solver.")
        @warn "Note: Lax-Friedrichs is diffusive and not able to capture shocks well..."
        println("Its stability is what makes it useful for checking results.")
    elseif solver_choice == 3
        println("You have chosen the Richtmyer solver.")
        @warn "Note: Richtmyer is a high-resolution method but can introduce spurious oscillations..."
        println("It is suitable for capturing shocks more accurately, with the oscillatiions as a caveat...")
    end
end

function ValidateGlobalUserInput(mode_choice, solver_choice, cfl, output_filename)
    output_filename = validate_mode_choice(mode_choice, output_filename)

    validate_solver_choice(solver_choice)

    if cfl <= 0
        println("CFL number must be positive. Defaulting to 0.5.")
        cfl = 0.5
    elseif cfl > 1
        @warn "CFL number should typically be less than or equal to 1 as it may lead to instability in some cases. Proceeding with the provided value."
    end
    return output_filename, cfl
end

function UI_sanity_check(UserInput)
    # Print selections
    println("\nYou selected:")
    for (k,v) in pairs(UserInput)
        println("  $k ---> $v")
    end
end