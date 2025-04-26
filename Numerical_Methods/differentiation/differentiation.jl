


include("derivatives_include.jl")

struct differentiation
    derivative::Symbol
    method::Symbol
    mode::Symbol
    accuracy::Union{Nothing, Int}
end

function setup_differentiation(derivative::Symbol, method::Symbol, mode::Symbol, accuracy::Int)
    return differntiation(derivative, method, mode, accuracy)
end


function forward_difference_dispatch(diff_config::differentiation, ind_var::AbstractVector{T}, dep_var::AbstractVector{T}, interp_ind_var) where T <: Real
    
    if diff_config.mode == :optimize

        coefficients = get_derivative_coefficients(diff_config.method, diff_config.derivative, diff_config.accuracy)

        if length(coefficients) != length(ind_var)

            interp_func(var) = interpolate_to_point(ind_var, dep_var, var)    
            interp_points = collect(range(ind_var[1], ind_var[end], length(coefficients)))
            h = (interp_points[end] - interp_points[1]) / (length(coefficients) - 1)
            dep_var = interp_func.(interp_points)  # Ensure broadcasting is applied properly

            return sum(coefficients .* y) /  h

        else

            h = (ind_var[end] - ind_var[1]) / (length(ind_var) - 1)
    
            h_optimal = sqrt(eps(Float64)) * abs(sum(ind_var)/length(ind_var))  

            if h < h_optimal

                @warn "Step size too small, potential numerical issues."
                h = h_optimal

            end

            return sum(coefficients .* y) /  h
            
        end   

    elseif diff_config.mode == :test
        # run_diff_convergence_test()
    elseif diff_config.mode == :debug
        # run_diff_debug_test()
    elseif diff_config.mode == :parallel
        # run_diff_parallel_test()
    else
        println("Unknown mode for forward difference")
    end 

    return forward_difference(ind_var, dep_var, interp_ind_var, diff_config.accuracy)
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
