




function charm_limiter(u, dx)
    # Number of grid points
    N = length(u)
    
    # Pre-allocate arrays for the limited solution and slopes
    u_limited = zeros(N)
    slopes = zeros(N)
    
    # Compute the slopes (differences) between adjacent points
    for i in 2:N-1
        slopes[i] = (u[i+1] - u[i-1]) / (2 * dx)
    end
    
    # Apply the CHARM limiter
    for i in 2:N-1
        # Calculate the limited slope
        r = (u[i] - u[i-1]) / (u[i+1] - u[i])
        phi = max(0, min(2 * r, 1))
        
        # Apply the limited slope to get the new value
        u_limited[i] = u[i] - 0.5 * phi * (u[i+1] - u[i-1])
    end
    
    return u_limited
end

# Example usage:
dx = 0.1
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # Sample solution vector

u_limited = charm_limiter(u, dx)
println("Limited Solution: ", u_limited)
