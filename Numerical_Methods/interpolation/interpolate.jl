


include("interpolation_include.jl")

struct interpolation
    method::Symbol
    mode::Symbol
    extrapolate::Bool
end

function setup_interpolation(method::Symbol, mode::Symbol, extrapolate::Bool)
    return interpolation(method, mode, extrapolate)
end


function interpolation_dispatch(ind_var::AbstractVector{T}, dep_var::AbstractVector{T}; interp_config::interpolation, d_dep_var=nothing) where T <: Real
    
    # Check if the lengths of ind_var and dep_var are equal
    if length(ind_var) != length(dep_var)
        throw(ArgumentError("ind_var and dep_var must be the same length"))
    end

    # Check if ind_var is sorted
    if !issorted(ind_var)
        throw(ArgumentError("ind_var must be sorted"))
    end

    # Check if ind_var and dep_var are not empty
    if isempty(ind_var) || isempty(dep_var)
        throw(ArgumentError("ind_var and dep_var must not be empty"))
    end
    
    # Check if ind_var and dep_var are of the same dimension
    if ndims(ind_var) != ndims(dep_var)
        throw(ArgumentError("ind_var and dep_var must be of the same dimension"))
    end

    # Check if ind_var and dep_var are of the same size
    if size(ind_var) != size(dep_var)
        throw(ArgumentError("ind_var and dep_var must be of the same size"))
    end

    if interp_config.method == :cubic_spline

        if interp_config.mode == :optimize
    
            if interp_config.extrapolate
                error("Extrapolation is not implemented yet")
            else
                return cubic_spline(ind_var, dep_var; throw_on_bounds = false)
            end
    
        elseif interp_config.mode == :test
            
            run_interp_convergence_test()
            run_interp_ptp_vs_all_points_test()

        elseif interp_config.mode == :debug
            
            return cubic_spline_debug(ind_var, dep_var; throw_on_bounds = false)

        elseif interp_config.mode == :parallel

            error("parallel mode not implemented yet")

        else
            println("Unknown mode for cubic spline interpolation")
        end
    
    elseif interp_config.method == :div_diff

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
            println("Unknown mode for cubic spline interpolation")
        end

        println("Using divided difference interpolation")
    elseif interp_config.method == :hermite

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
            
            #run_interp_convergence_test()
            #run_interp_ptp_vs_all_points_test()
        elseif interp_config.mode == :debug
            
            return hermite_interpolation_logged(ind_var, dep_var, d_dep_var, point)

        elseif interp_config.mode == :parallel

    elseif interp_config.method == :WENO
        
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

        else
            println("Unknown method for interpolation")
        end
    else
        println("Unknown interpolation method")
    end
end     




# Choices for method are :cubic_spline, :div_diff, :hermite, :WENO 
method = :cubic_spline
# Choices for mode are :optimize, :test, :debug, :parallel 
mode = :optimize
# Choices for eind_vartrapolation are true or false
eind_vartrapolate = false
interpolate = setup_interpolation(method, mode, extrapolate)

ind_var = collect(range(0.0, stop=10.0, length=100))
dep_var = sin.(ind_var)



println("Interpolation method: ", interpolate.method)
println("Interpolation mode: ", interpolate.mode)
println("Extrapolation: ", interpolate.extrapolate)