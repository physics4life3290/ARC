using Plots

"""
Minerbo Closure (Variable Eddington Factor)
"""
function minerbo_closure(E, F)
    f = abs(E) > 1e-10 ? min(abs(F) / E, 0.9999) : 0.0
    return (1/3) + (2 * f^2) / (3 * (1 + sqrt(1 - f^2)))
end

"""
HLLE Flux Solver
Ensures the flux between cells is physically consistent and stable.
"""
function hlle_flux(EL, FL, ER, FR)
    c = 1.0
    # Characteristic speeds (approximate for M1)
    lambda_L = -c
    lambda_R = c
    
    # Pressure from closure
    PL = minerbo_closure(EL, FL) * EL
    PR = minerbo_closure(ER, FR) * ER
    
    # State vectors and Flux vectors
    # U = [E, F], Flux(U) = [F, c^2 * P]
    flux_E_L, flux_F_L = FL, c^2 * PL
    flux_E_R, flux_F_R = FR, c^2 * PR
    
    # HLLE Formula: F = (lambda_R*F_L - lambda_L*F_R + lambda_R*lambda_L*(U_R - U_L)) / (lambda_R - lambda_L)
    fE = (lambda_R * flux_E_L - lambda_L * flux_E_R + lambda_R * lambda_L * (ER - EL)) / (lambda_R - lambda_L)
    fF = (lambda_R * flux_F_L - lambda_L * flux_F_R + lambda_R * lambda_L * (FR - FL)) / (lambda_R - lambda_L)
    
    return fE, fF
end

function simulate_stable_transport()
    # Grid
    nx = 200
    dx = 0.05
    dt = 0.01 
    t_steps = 500
    c = 1.0

    # Fields
    E = zeros(nx)
    F = zeros(nx)
    
    # Initial Condition: A pulse of neutrinos
    x_axis = range(0, stop=(nx-1)*dx, length=nx)
    E .= [exp(-0.5 * ((x - 2.0) / 0.5)^2) for x in x_axis]
    
    # Opacity (Matter Background)
    kappa = zeros(nx)
    kappa[120:160] .= 10.0 # A "wall" of matter

    for t in 1:t_steps
        E_new = copy(E)
        F_new = copy(F)
        
        # We calculate fluxes at the interfaces (i + 1/2)
        # Using a loop to update interior cells
        for i in 2:nx-1
            # HLLE Fluxes at left and right interfaces
            fE_left, fF_left = hlle_flux(E[i-1], F[i-1], E[i], F[i])
            fE_right, fF_right = hlle_flux(E[i], F[i], E[i+1], F[i+1])
            
            # Source terms (Absorption)
            S_E = -kappa[i] * c * E[i]
            S_F = -kappa[i] * c * F[i]
            
            # Finite Volume Update
            E_new[i] += dt * (-(fE_right - fE_left) / dx + S_E)
            F_new[i] += dt * (-(fF_right - fF_left) / dx + S_F)
        end
        
        E, F = E_new, F_new
        
        # Simple vacuum boundary conditions (Outflow)
        E[1] = E[2]; E[nx] = E[nx-1]
        F[1] = F[2]; F[nx] = F[nx-1]
    end
    
    return x_axis, E, kappa
end

# Execution
x, E_final, k = simulate_stable_transport()
p1 = plot(x, E_final, label="Neutrino Energy", color=:blue, lw=2)
p1 = plot!(x, k ./ maximum(k), label="Matter Opacity (scaled)", color=:gray, fill=(0, 0.2, :gray))
display(p1)
