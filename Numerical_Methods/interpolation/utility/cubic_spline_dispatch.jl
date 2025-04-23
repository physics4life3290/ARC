




function cubic_spline_dispatch(interp_config, ind_var, dep_var, interp_ind_var)
    
    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return cubic_spline(ind_var, dep_var, interp_ind_var; throw_on_bounds = false)
        end

    elseif interp_config.mode == :test
        
        interp_fn(x, y, interp_ind_var) = [cubic_spline(x, y, xi; throw_on_bounds=false) for xi in interp_ind_var]
        run_interp_convergence_test(interp_fn)

        #run_interp_ptp_vs_all_points_test()

    elseif interp_config.mode == :debug
        
        return cubic_spline_debug(ind_var, dep_var, interp_ind_var)

    elseif interp_config.mode == :parallel

        error("parallel mode not implemented yet")

    else
        println("Unknown mode for cubic spline interpolation")
    end
end