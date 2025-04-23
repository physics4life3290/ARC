




function hermite_dispatch(interp_config, ind_var, dep_var, d_dep_var, point)     

    if d_dep_var === nothing
        throw(ArgumentError("d_dep_var must be provided for Hermite interpolation"))
    end

    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return hermite_interpolation(ind_var, dep_var, d_dep_var, point)
        end
    elseif interp_config.mode == :test
        
        interp_fn(ind_var, dep_var, d_dep_var, interp_ind_var) = [hermite_interpolation(ind_var, dep_var, d_dep_var, xi) for xi in interp_ind_var]
        run_interp_convergence_test(interp_fn)
        #run_interp_ptp_vs_all_points_test()
        
    elseif interp_config.mode == :debug
        
        return hermite_interpolation_logged(ind_var, dep_var, d_dep_var, point)

    elseif interp_config.mode == :parallel
            
            error("parallel mode not implemented yet")
    else
        println("Unknown mode for Hermite interpolation")
    end
end