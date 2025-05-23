




function div_diff_dispatch(interp_config, ind_var, dep_var, interp_ind_var)

    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return newton_interpolation(ind_var, dep_var, interp_ind_var)
        end

    elseif interp_config.mode == :test
        
        run_interp_convergence_test(newton_interpolation)
        
    elseif interp_config.mode == :debug
        
        return newton_interpolation_logged(ind_var, dep_var, interp_ind_var)

    elseif interp_config.mode == :parallel

        error("parallel mode not implemented yet")

    else
        error("Unknown mode for Divided Difference interpolation")
    end
end