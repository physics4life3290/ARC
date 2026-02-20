




# Compute slopes for primitive variables with minmod limiter
function compute_slopes!(slopes, var, limiter)
    n = length(var)

    slopes[1] = 0.0
    slopes[n] = 0.0

    @inbounds for i in 2:n-1
        dl = var[i] - var[i-1]
        dr = var[i+1] - var[i]

        if limiter === :minmod
            slopes[i] = minmod(dl, dr)
        elseif limiter === :superbee
            slopes[i] = superbee(dl, dr)
        elseif limiter === :vanleer
            slopes[i] = vanleer(dl, dr)
        else
            error("Unknown limiter")
        end
    end

    return nothing
end

function linear_reconstruct(var, slopes)
    n = length(var)

    var_L = zeros(n + 1)
    var_R = zeros(n + 1)

    # interior interfaces
    @inbounds for i in 1:n-1
        var_L[i+1] = var[i]   + 0.5 * slopes[i]
        var_R[i+1] = var[i+1] - 0.5 * slopes[i+1]
    end

    # boundaries: piecewise constant
    var_L[1]   = var[1]
    var_R[1]   = var[1]

    var_L[n+1] = var[n]
    var_R[n+1] = var[n]

    return var_L, var_R
end
