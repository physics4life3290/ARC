




function cubic_spline_dispatch(interp_config, ind_var, dep_var, interp_ind_var)
    
    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return cubic_spline(ind_var, dep_var, interp_ind_var; throw_on_bounds = false)
        end

    elseif interp_config.mode == :test
        
        interpolation_convergence_test(cubic_spline)

    elseif interp_config.mode == :debug
        
        return cubic_spline_debug(ind_var, dep_var, interp_ind_var)

    elseif interp_config.mode == :parallel

        error("parallel mode not implemented yet")

    else
        println("Unknown mode for cubic spline interpolation")
    end
end