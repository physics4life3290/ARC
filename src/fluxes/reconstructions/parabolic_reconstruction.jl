




function parabolic_reconstruct(var::AbstractVector)
    n = length(var)
    var_L = zeros(n)
    var_R = zeros(n)

    # Interior cells
    Threads.@threads for i in 2:(n-1)
        @inbounds begin
            d_im1 = var[i] - var[i-1]
            d_ip1 = var[i+1] - var[i]

            slope = 0.5 * (d_im1 + d_ip1)          # first derivative
            curvature = 0.5 * (d_ip1 - d_im1)     # second derivative

            # Left and right interface reconstruction
            var_R[i] = var[i] + 0.5*slope + (1/6)*curvature   # right interface i+1/2
            var_L[i] = var[i] - 0.5*slope + (1/6)*curvature   # left interface i-1/2
        end
    end

    # Boundaries: copy cell-centered values
    var_L[1] = var[1]; var_R[1] = var[1]
    var_L[end] = var[end]; var_R[end] = var[end]

    return var_L, var_R
end