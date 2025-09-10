




function parabolic_reconstruct(var::AbstractArray)
    n = length(var)
    # Allocate left/right interface arrays
    var_L = zeros(n)
    var_R = zeros(n)

    # Parabolic reconstruction coefficients:
    # c0 = cell average (var[i])
    # c1 = slope term
    # c2 = curvature term

    Threads.@threads for i in 2:(n-1)
        @inbounds begin
            # Compute first differences
            d_im1 = var[i] - var[i-1]
            d_ip1 = var[i+1] - var[i]

            # Centered slope (first derivative)
            slope = 0.5 * (d_im1 + d_ip1)

            # Quadratic curvature (second derivative)
            curvature = 0.5 * (d_ip1 - d_im1)

            # Reconstruct left and right interface states
            # L state at i+1/2 comes from cell i
            var_R[i] = var[i] + 0.5 * slope + (1/6) * curvature

            # R state at i-1/2 comes from cell i
            var_L[i] = var[i] - 0.5 * slope + (1/6) * curvature
        end
    end

    # Boundaries: copy cell values
    var_L[1] = var[1]
    var_R[1] = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end