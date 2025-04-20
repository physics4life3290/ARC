# Define the RK4 solver function
function rk4(f, u0, tspan, dt)
    t0, tf = tspan
    t = t0:dt:tf  # Time vector from t0 to tf with step size dt
    u = zeros(length(t), length(u0))  # Array to store the solutions
    u[1, :] = u0  # Initial condition
    
    # Iterate over time steps
    for i in 1:(length(t) - 1)
        k1 = dt * f(t[i], u[i, :])
        k2 = dt * f(t[i] + dt/2, u[i, :] + k1/2)
        k3 = dt * f(t[i] + dt/2, u[i, :] + k2/2)
        k4 = dt * f(t[i] + dt, u[i, :] + k3)
        
        # Update the solution
        u[i + 1, :] = u[i, :] + (k1 + 2*k2 + 2*k3 + k4) / 6
    end
    
    return t, u
end

# Example usage: Solving the simple ODE du/dt = -u with initial condition u(0) = 1
f(t, u) = -ℯ.^-u .* cos.(u)  # Define the differential equation
u0 = [1.0]  # Initial condition
tspan = (0.0, 10.0)  # Time range
dt = 0.1  # Time step size

# Solve the ODE
t, u = rk4(f, u0, tspan, dt)

# Plot the results (optional)
using Plots
plot(t, u[:, 1], label="u(t)", xlabel="t", ylabel="u")
