




function weno_dispatch(interp_config, ind_var::AbstractVector{T}, dep_var::AbstractVector{T}; point::Real) where T <: Real
    if interp_config.mode == :optimize

        if interp_config.extrapolate
            error("Extrapolation is not implemented yet")
        else
            return WENO(ind_var, dep_var)
        end

    elseif interp_config.mode == :test
        
        #run_interp_convergence_test()
        #run_interp_ptp_vs_all_points_test()

    elseif interp_config.mode == :debug
        
        return  weno5_interpolate_at_logged(ind_var, dep_var, point; interpolation_type = :left, logfile = "Numerical_Methods/interpolation/debug/output/weno5_debug.txt")

    elseif interp_config.mode == :parallel

        error("parallel mode not implemented yet")

    else
        println("Unknown mode for WENO interpolation")
    end
end