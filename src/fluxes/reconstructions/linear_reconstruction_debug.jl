


function compute_slopes_debug(var::AbstractArray, n::Int; limiter_input::Symbol=:minmod, mode=:serial)

    exceptions_serial = Atomic{Union{Nothing, Exception}}(nothing)
    exceptions_threads = Threads.Atomic{Union{Nothing, Exception}}(nothing)

    slopes = zeros(n)

    if mode == :serial

        for i in 2:(n-1)
            try
                @inbounds begin
                    dl = var[i] - var[i-1]
                    dr = var[i+1] - var[i]
                    if limiter_input !== nothing 
                        if limiter_input == :minmod
                            slopes[i] = minmod(dl, dr)
                        elseif limiter_input == :superbee
                            slopes[i] = superbee(dl, dr)
                        elseif limiter_input == :vanleer
                            slopes[i] = vanleer(dl, dr)
                        end
                    end 
                end
            catch e
                atomic_compare_exchange!(exceptions_serial, nothing, e)
            end
        end

    elseif mode == :threads 

        Threads.@threads for i in 2:(n-1)
            try
                @inbounds begin
                    dl = var[i] - var[i-1]
                    dr = var[i+1] - var[i]
                    if limiter_input !== nothing 
                        if limiter_input == :minmod
                            slopes[i] = minmod(dl, dr)
                        elseif limiter_input == :superbee
                            slopes[i] = superbee(dl, dr)
                        elseif limiter_input == :vanleer
                            slopes[i] = vanleer(dl, dr)
                        end
                    end
                end
            catch e 
                atomic_compare_exchange!(exceptions_threads, nothing, e)
            end
        end

    end

    slopes[1] = 0.0
    slopes[end] = 0.0

    return slopes
end

# Reconstruction at interfaces: left and right states
function linear_reconstruct(var::AbstractArray, slopes::AbstractArray, n::Int; mode=:serial)

    exceptions_serial = Atomic{Union{Nothing, Exception}}(nothing)
    exceptions_threads = Threads.Atomic{Union{Nothing, Exception}}(nothing)

    # Left state at interface i+1/2 comes from cell i
    var_L = zeros(n)
    # Right state at interface i+1/2 comes from cell i+1
    var_R = zeros(n)

    if mode == :serial

        for i in 1:n
            try
                @inbounds begin
                    var_L[i] = var[i] + 0.5 * slopes[i]
                    var_R[i] = var[i] - 0.5 * slopes[i]
                end
            catch e 
                atomic_compare_exchange!(exceptions_serial, nothing, e)
            end
        end

    elseif mode == :thread

        # Interior interfaces
        Threads.@threads for i in 1:n
            try
                @inbounds begin
                    var_L[i] = var[i] + 0.5 * slopes[i]
                    var_R[i] = var[i] - 0.5 * slopes[i]
                end
            catch e 
                atomic_compare_exchange!(exceptions_threads, nothing, e)
            end
        end
        
    end

    return var_L, var_R
end