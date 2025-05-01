


include("interpolation_include.jl")

struct interpolation
    method::Symbol
    mode::Symbol
    extrapolate::Bool
end

function setup_interpolation(method::Symbol, mode::Symbol, extrapolate::Bool)
    return interpolation(method, mode, extrapolate)
end


function interpolation_dispatch(ind_var::AbstractVector{T}, dep_var::AbstractVector{T}, interp_ind_var; interp_config=:default, d_dep_var=nothing) where T <: Real
    if interp_config == :default
        method = :cubic_spline
        mode = :optimize
        extrapolate = false
        interp_config = setup_interpolation(method, mode, extrapolate)
    end
    # Check input for proper characteristics
    prelim_interp_check(ind_var, dep_var)
 
    if interp_config.method == :cubic_spline

        return cubic_spline_dispatch(interp_config, ind_var, dep_var, interp_ind_var)
    
    elseif interp_config.method == :div_diff
        @warn("Warning: This method works best when used piecewise. This is due to Gibb's Phenomena!")
        return div_diff_dispatch(interp_config, ind_var, dep_var, interp_ind_var)
        
    elseif interp_config.method == :hermite
        
        @warn("Warning: This method works best when used piecewise. This is due to Gibb's Phenomena!")
        return hermite_dispatch(interp_config, ind_var, dep_var, d_dep_var, interp_ind_var)

    elseif interp_config.method == :WENO
        
        return weno_dispatch(interp_config, ind_var, dep_var, interp_ind_var)

    else
        println("Unknown interpolation method")
    end
end     



#=
# Choices for method are :cubic_spline, :div_diff, :hermite, :WENO 
method = :div_diff
# Choices for mode are :optimize, :test, :debug, :parallel 
mode = :test
# Choices for extrapolation are true or false
extrapolate = false
interpolate = setup_interpolation(method, mode, extrapolate)
ind_var = collect(range(0.0, stop=10.0, length=100))
#dep_var = sin.(ind_var)
#d_dep_var = cos.(ind_var)
dep_var = ind_var .^ 2
d_dep_var = 2 .* ind_var
point = 1.75

println(interpolation_dispatch(ind_var, dep_var, point; interp_config=interpolate, d_dep_var=d_dep_var))
=#