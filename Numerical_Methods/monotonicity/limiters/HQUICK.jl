




# HQUICK Limiter function
function hquick_limiter(u, dx)
    # Number of grid points
    N = length(u)
    
    # Pre-allocate arrays for the limited solution and slopes
    u_limited = zeros(N)
    
    # Function for QUICK interpolation (quadratic interpolation)
    function quick_interpolation(i)
        return 1.5 * u[i] - 0.5 * u[i+1] + 0.5 * (u[i+1] - u[i-1])
    end
    
    # Apply the HQUICK limiter
    for i in 2:N-1
        # Compute the unconstrained QUICK interpolation
        q = quick_interpolation(i)
        
        # Compute the local slopes
        slope1 = (u[i] - u[i-1]) / dx
        slope2 = (u[i+1] - u[i]) / dx
        
        # Determine the limiting factor (phi) based on the ratio of slopes
        r = min(slope1, slope2)
        
        # Use the ratio r to constrain the interpolation if necessary
        phi = max(0, min(2 * r, 1))  # Constrained by a limiting factor
        
        # Apply the limited interpolation to get the new value
        u_limited[i] = u[i] + phi * (q - u[i])
    end
    
    return u_limited
end

# Example usage:
dx = 0.1
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # Sample solution vector

# Apply HQUICK limiter
u_limited = hquick_limiter(u, dx)
println("Limited Solution using HQUICK: ", u_limited)
