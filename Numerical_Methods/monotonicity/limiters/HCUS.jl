# HCUS Limiter function with method as a positional argument
function hcus_limiter(u, dx, method)
    # Number of grid points
    N = length(u)
    
    # Pre-allocate arrays for the limited solution and slopes
    u_limited = zeros(N)
    slopes = zeros(N)
    
    # Compute the slopes (differences) between adjacent points
    for i in 2:N-1
        slopes[i] = (u[i+1] - u[i-1]) / (2 * dx)
    end
    
    # Function for the minmod limiter (one of the most common choices)
    function minmod(a, b)
        return sign(a) * min(abs(a), abs(b))
    end
    
    # Function for the Superbee limiter
    function superbee(r)
        return max(0.0, min(2 * r, 1.0), min(r, 2.0))
    end
    
    # Apply the HCUS limiter using the selected method
    for i in 2:N-1
        # Calculate the local ratio for the constrained slope limiter
        r = (u[i] - u[i-1]) / (u[i+1] - u[i])
        
        # Select the limiting function
        if method == :minmod
            phi = minmod(slopes[i], slopes[i-1])
        elseif method == :superbee
            phi = superbee(r)
        else
            error("Unsupported limiter method")
        end
        
        # Apply the limited slope to get the new value
        u_limited[i] = u[i] - 0.5 * phi * (u[i+1] - u[i-1])
    end
    
    return u_limited
end

# Example usage:
dx = 0.1
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # Sample solution vector

# Apply HCUS limiter using Minmod method
u_limited_minmod = hcus_limiter(u, dx, :minmod)
println("Limited Solution using Minmod: ", u_limited_minmod)

# Apply HCUS limiter using Superbee method
u_limited_superbee = hcus_limiter(u, dx, :superbee)
println("Limited Solution using Superbee: ", u_limited_superbee)