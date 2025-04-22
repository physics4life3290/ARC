




function div_diff_dispatch(interp_config, ind_var, dep_var)

    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return div_diff(ind_var, dep_var)
        end

    elseif interp_config.mode == :test
        
        #run_interp_convergence_test()
        #run_interp_ptp_vs_all_points_test()

    elseif interp_config.mode == :debug
        
        return newton_divided_differences_logged(ind_var, dep_var)

    elseif interp_config.mode == :parallel

        error("parallel mode not implemented yet")

    else
        error("Unknown mode for Divided Difference interpolation")
    end
end