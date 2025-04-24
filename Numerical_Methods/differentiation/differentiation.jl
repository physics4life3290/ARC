




struct differentiation
    derivative::Symbol
    method::Symbol
    mode::Symbol
    accuracy::Int
    interpolate::Bool
end

function setup_differentiation(derivative::Symbol, method::Symbol, mode::Symbol, accuracy::Int, interpolate::Bool)
    return differntiation(derivative, method, mode, accuracy, interpolate)
end

function differentiation_dispatch(ind_var::AbstractVector{T}, dep_var::AbstractVector{T}, interp_ind_var; diff_config::differentiation, d_dep_var=nothing) where T <: Real
    
    # Check input for proper characteristics
    prelim_diff_check(ind_var, dep_var)

    if diff_config.derivative == :forward
        return forward_difference_dispatch(diff_config, ind_var, dep_var, interp_ind_var)
    elseif diff_config.derivative == :backward
        return backward_difference_dispatch(diff_config, ind_var, dep_var, interp_ind_var)
    elseif diff_config.derivative == :central
        return central_difference_dispatch(diff_config, ind_var, dep_var, interp_ind_var)
    else
        println("Unknown derivative method")
    end
end
