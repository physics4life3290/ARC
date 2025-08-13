




# Compute slopes for primitive variables with minmod limiter
function compute_slopes(var, limiter_input)
    n = length(var)
    slopes = zeros(n)
    Threads.@threads for i in 2:n-1
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
    end
    slopes[1] = 0.0
    slopes[end] = 0.0
    return slopes
end

# Reconstruction at interfaces: left and right states
function linear_reconstruct(var, slopes)

    n = length(var)
    # Left state at interface i+1/2 comes from cell i
    var_L = zeros(n)
    # Right state at interface i+1/2 comes from cell i+1
    var_R = zeros(n)

    # Interior interfaces
    Threads.@threads for i in 2:n
        @inbounds begin
            var_L[i] = var[i] + 0.5 * slopes[i]
            var_R[i] = var[i] - 0.5 * slopes[i]
        end
    end

    # Boundary interfaces (use cell values, no extrapolation)
    var_L[1] = var[1]
    var_R[1] = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end