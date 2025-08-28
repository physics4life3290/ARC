




function flattening_curvature_debug(ρ::AbstractVector, slope_ρ::AbstractVector, n::Int, dx::Real)

    exceptions_serial = Atomic{Union{Nothing, Exception}}(nothing)
    exceptions_threads = Threads.Atomic{Union{Nothing, Exception}}(nothing)
    # We should say something in debug about only density quantities are used here.

    for i in 2:(n-1)

        try
            @inbounds begin
                # First derivative (centered)
                f_prime = (ρ[i+1] - ρ[i-1]) / (2*dx)

                # Second derivative (centered)
                f_double_prime = (ρ[i+1] - 2ρ[i] + ρ[i-1]) / (dx^2)

                # True curvature
                curvature = abs(f_double_prime) / (1 + f_prime^2)^(3/2)

                # Flattening coefficient
                flattening_coef = clamp(1.0 - 10.0 * curvature, 0.0, 1.0)

                slope_ρ[i] *= flattening_coef
            end
        catch e
            atomic_compare_exchange!(exceptions_serial, nothing, e)
        end

        if verbose == true
            @info "i=$i, f'=$f_prime, f''=$f_double_prime, curvature=$curvature, flattening_coef=$flattening_coef"
        end

    end

    Threads.@threads for i in 2:(n-1)

        try
            @inbounds begin
                f_prime = (ρ[i+1] - ρ[i-1]) / (2*dx)
                f_double_prime = (ρ[i+1] - 2*ρ[i] + ρ[i-1]) / (dx^2)
                curvature = abs(f_double_prime) / (1 + f_prime^2)^(3/2)
                flattening_coef = clamp(1.0 - 10.0 * curvature, 0.0, 1.0)
                slope_ρ[i] *= flattening_coef
            end
        catch e
            # Store the first exception encountered
            Threads.atomic_compare_exchange!(exceptions, nothing, e)
        end

        if verbose == true
            @info "i=$i, f'=$f_prime, f''=$f_double_prime, curvature=$curvature, flattening_coef=$flattening_coef"
        end

    end

    # If an exception occurred, throw it
    if exceptions_serial[] !== nothing
        throw(exceptions_serial[])
    elseif exceptions_threads[] !== nothing
        throw(exceptions_threads[])
    end

    return slope_ρ
end



function compute_flattening_coefficient_debug(var::AbstractVector, n::Int; ε=1e-10)

    θ = ones(n)  # start with no flattening (1.0)
    exceptions_serial = Atomic{Union{Nothing, Exception}}(nothing)
    exceptions_threads = Threads.Atomic{Union{Nothing, Exception}}(nothing)

    Threads.@threads for i in 2:(n-1)
        try
            @inbounds begin
                Δp_plus = abs(var[i+1] - var[i])
                Δp_minus = abs(var[i] - var[i-1])
                Δp = max(Δp_plus, Δp_minus)
                p_avg = (abs(var[i+1]) + 2*abs(var[i]) + abs(var[i-1])) / 4

                # Avoid division by zero
                shock_strength = Δp / (p_avg + ε)

                # Thresholds for flattening
                if shock_strength > 0.1
                    θ[i] = 0.0  # full flattening near strong shock
                elseif shock_strength > 0.05
                    θ[i] = 0.5  # partial flattening near moderate gradient
                else
                    θ[i] = 1.0  # no flattening in smooth regions
                end
            end
        catch e
            Threads.atomic_compare_exchange!(exceptions, nothing, e)
        end
    end

    # If an exception occurred, throw it
    if exceptions_serial[] !== nothing
        throw(exceptions_serial[])
    elseif exceptions_threads[] !== nothing
        throw(exceptions_threads[])
    end

    # Boundaries: replicate neighbors
    θ[1] = θ[2]
    θ[end] = θ[end-1]

    return θ
end

function compute_flattening_coefficient_debug(var::AbstractVector, n::Int; ε=1e-10; mode::Symbol=:serial, verbose::Bool=false)
    
    θ = ones(n)  # start with no flattening (1.0)

    if mode == :serial

        for i in 2:(n-1)
            @inbounds begin
                Δp_plus = abs(var[i+1] - var[i])
                Δp_minus = abs(var[i] - var[i-1])
                Δp = max(Δp_plus, Δp_minus)
                p_avg = (abs(var[i+1]) + 2*abs(var[i]) + abs(var[i-1])) / 4

                # Avoid division by zero, ε small number
                shock_strength = Δp / (p_avg + ε)

                # Thresholds for flattening (tweak these!)
                if shock_strength > 0.1
                    θ[i] = 0.0  # full flattening near strong shock
                elseif shock_strength > 0.05
                    θ[i] = 0.5  # partial flattening near moderate gradient
                else
                    θ[i] = 1.0  # no flattening in smooth regions
                end
                if verbose == true
                    @info "i=$i, Δp_plus=$Δp_plus, Δp_minus=$Δp_minus, Δp=$Δp, p_avg=$p_avg, shock_strength=$shock_strength, θ[i]=${θ[i]}"
                end
            end
        end 

    elseif mode == :threads

        Threads.@threads for i in 2:(n-1)
            @inbounds begin
                Δp_plus = abs(var[i+1] - var[i])
                Δp_minus = abs(var[i] - var[i-1])
                Δp = max(Δp_plus, Δp_minus)
                p_avg = (abs(var[i+1]) + 2*abs(var[i]) + abs(var[i-1])) / 4

                # Avoid division by zero, ε small number
                shock_strength = Δp / (p_avg + ε)

                # Thresholds for flattening (tweak these!)
                if shock_strength > 0.1
                    θ[i] = 0.0  # full flattening near strong shock
                elseif shock_strength > 0.05
                    θ[i] = 0.5  # partial flattening near moderate gradient
                else
                    θ[i] = 1.0  # no flattening in smooth regions
                end
                if verbose == true
                    @info "i=$i, Δp_plus=$Δp_plus, Δp_minus=$Δp_minus, Δp=$Δp, p_avg=$p_avg, shock_strength=$shock_strength, θ[i]=$(θ[i])"
                end
            end
        end

    end
    # Boundaries: no flattening (or replicate neighbors)
    θ[1] = θ[2]
    θ[end] = θ[end-1]

    return θ
  

end