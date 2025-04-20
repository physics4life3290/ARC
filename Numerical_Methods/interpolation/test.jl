


include("interpolation_include.jl")

struct interpolation
    method::Symbol
    mode::Symbol
    eind_vartrapolate::Bool
end

function setup_interpolation(method::Symbol, mode::Symbol, eind_vartrapolate::Bool)
    return interpolation(method, mode, eind_vartrapolate)
end


function interpolation_dispatch(interp_config::interpolation, ind_var::AbstractVector{T}, dep_var::AbstractVector{T}) where T <: Real
    
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
            
            error("Testing not implemented yet")

        elseif interp_config.mode == :debug
            
            error("debug mode not implemented yet")

        elseif interp_config.mode == :parallel

            error("parallel mode not implemented yet")

        else
            println("Unknown mode for cubic spline interpolation")
        end
    
    elseif interp_config.method == :div_diff
        println("Using divided difference interpolation")
    elseif interp_config.method == :hermite
        println("Using Hermite interpolation")
    elseif interp_config.method == :WENO
        println("Using WENO interpolation")
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
interpolate = setup_interpolation(method, mode, eind_vartrapolate)

ind_var = collect(range(0.0, stop=10.0, length=100))
dep_var = sin.(ind_var)



println("Interpolation method: ", interpolate.method)
println("Interpolation mode: ", interpolate.mode)
println("Eind_vartrapolation: ", interpolate.eind_vartrapolate)