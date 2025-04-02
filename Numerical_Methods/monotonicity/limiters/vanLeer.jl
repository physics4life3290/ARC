




# Van Leer Limiter function
function van_leer_limiter(u, dx)
    # Number of grid points
    N = length(u)
    
    # Pre-allocate array for the limited solution
    u_limited = zeros(N)
    
    # Function to calculate the ratio of slopes
    function calc_ratio(i)
        return (u[i+1] - u[i]) / (u[i] - u[i-1])
    end
    
    # Van Leer limiter function
    function van_leer_phi(r)
        return (r + abs(r)) / (1 + abs(r))
    end
    
    # Apply the Van Leer limiter
    for i in 2:N-1
        # Calculate the slope ratio
        r = calc_ratio(i)
        
        # Calculate the limiting factor using Van Leer's function
        phi = van_leer_phi(r)
        
        # Apply the limiter to the solution
        u_limited[i] = u[i] - 0.5 * phi * (u[i+1] - u[i-1])
    end
    
    return u_limited
end

# Example usage:
dx = 0.1
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # Sample solution vector

# Apply Van Leer limiter
u_limited = van_leer_limiter(u, dx)
println("Limited Solution using Van Leer: ", u_limited)
