




function collect_solver_input(primary_input)

    limiters = [:Minmod, :Superbee, :VanLeer]
    reconstructions = [:Constant, :Linear, :Parabolic]
    riemann_solver_choices = [:HLL, :HLLC, :Exact]
    bool_choices = [true, false]
    
    if primary_input.dimension == 1
        println("Please select desired options for $(primary_input.solver)...")
        cfl = prompt_cfl()
        t_final = prompt_t_final()
        
        if primary_input.solver == :FTCS || primary_input.solver == :LaxFriedrichs || primary_input.solver == :Richtmyer
            return SolverInput(cfl, t_final, false, false, nothing, nothing, nothing)
        elseif primary_input.solver == :GodunovScheme
            println("Which interface reconstruction?")
            reconstruction = prompt_choice("Choose a reconstruction method: ", reconstructions)
            if reconstruction == :Constant
                @warn "Flattening, Steepening, and Limiters are not applicable for Constant Reconstruction and will be ignored..."
            elseif reconstruction == :Linear || reconstruction == :Parabolic
                println("Would you like flattening of the shock waves?")
                flattening = prompt_choice_bool("Choose true of false: ", bool_choices)
                println("Would you like steepening of the contact discontinuities?")
                steepening = prompt_choice_bool("Choose true or false: ", bool_choices)
                println("Would you like to use a limiter?")
                limiter_choice = prompt_choice_bool("Choose true or false: ", bool_choices)
                if limiter_choice == true
                    limiter = prompt_choice("Choose a limiter:", limiters)
                else
                    limiter = nothing
                end
            end
            println("Which Riemann Solver would you like?")
            riemann_solve = prompt_choice("Choose a Riemann Solver: ", riemann_solver_choices)
            return SolverInput(cfl, t_final, flattening, steepening, limiter, reconstruction, riemann_solve)
        end
    end
end