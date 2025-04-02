




# OSPRE Limiter function
function ospre_limiter(u, dx)
    # Number of grid points
    N = length(u)
    
    # Pre-allocate arrays for the limited solution and slopes
    u_limited = zeros(N)
    
    # Function for calculating the slope between points
    function calc_slope(i)
        return (u[i+1] - u[i]) / dx
    end
    
    # Function to calculate the "smoothness" factor
    function smoothness_factor(i)
        return abs((u[i+1] - u[i]) / (u[i] - u[i-1]))
    end
    
    # Apply the OSPRE limiter
    for i in 2:N-1
        # Calculate the slope at each point
        slope = calc_slope(i)
        
        # Calculate the smoothness factor for the interpolation
        smoothness = smoothness_factor(i)
        
        # Limiting based on smoothness and slope
        if smoothness < 1.0  # If the solution is smooth
            u_limited[i] = u[i] + slope * (u[i+1] - u[i-1]) * 0.5  # Apply a conservative adjustment
        else  # If the solution has oscillations
            u_limited[i] = u[i]  # Keep the value as it is (no adjustment)
        end
    end
    
    return u_limited
end

# Example usage:
dx = 0.1
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # Sample solution vector

# Apply OSPRE limiter
u_limited = ospre_limiter(u, dx)
println("Limited Solution using OSPRE: ", u_limited)
