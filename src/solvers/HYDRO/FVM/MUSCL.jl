




function reconstruct_parabolic(var)
    n = length(var)

    # Left and Right reconstructed states at interfaces (n+1 interfaces)
    var_L = zeros(n+1)
    var_R = zeros(n+1)

    # Compute parabolic interface values using a 3-point stencil
    for i in 2:n-1
        # Fit a parabola through var[i-1], var[i], var[i+1]
        dq_minus = var[i] - var[i-1]
        dq_plus = var[i+1] - var[i]

        # Centered difference with limiter (basic minmod to reduce overshoot)
        slope = minmod(dq_minus, dq_plus)

        # Compute curvature
        d2q = var[i+1] - 2*var[i] + var[i-1]

        # Parabolic reconstruction
        var_L[i+1] = var[i] + 0.5 * slope - (1/6) * d2q
        var_R[i]   = var[i] - 0.5 * slope - (1/6) * d2q
    end

    # Boundary extrapolation (copy cell-centered value)
    var_L[1]   = var[1]
    var_R[1]   = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end

function GodunovStep!(user_input, _grid, W, U, dt)
    
    nx = _grid.xcoord.total_zones
    γ = user_input.Secondary_Input.gamma
    dx = _grid.xcoord.spacing

    if user_input.Solver_Input.reconstruction != :Constant
        # Compute slopes for primitive variables
        slope_ρ = compute_slopes(W.density_centers, user_input.Solver_Input.limiter)
        slope_u = compute_slopes(W.velocity_centers, user_input.Solver_Input.limiter)
        slope_p = compute_slopes(W.pressure_centers, user_input.Solver_Input.limiter)
        
        if user_input.Solver_Input.flattening == true
            flattening(W, slope_ρ)
        end

        if user_input.Solver_Input.steepening == true
            slope_ρ = steepening(W.density_centers, W.pressure_centers, slope_ρ, _grid.xcoord.total_zones)
        end
        
        # Reconstruct left and right states at interfaces
        if user_input.Solver_Input.reconstruction == :Linear
            ρL, ρR = linear_reconstruct(W.density_centers, slope_ρ)
            uL, uR = linear_reconstruct(W.velocity_centers, slope_u)
            pL, pR = linear_reconstruct(W.pressure_centers, slope_p)
        end

    else
        ρL, ρR = W.density_centers, W.density_centers
        uL, uR = W.velocity_centers, W.velocity_centers
        pL, pR = W.pressure_centers, W.pressure_centers
    end

    apply_boundary_conditions(user_input, U, _grid)
    # Fluxes at interfaces
    F1 = zeros(nx+1)
    F2 = zeros(nx+1)
    F3 = zeros(nx+1)
    Threads.@threads for i in 2:nx
        @inbounds begin
            if user_input.Solver_Input.riemanntype == :HLL
                f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i])
            elseif user_input.Solver_Input.riemanntype == :HLLC
                f1, f2, f3 = Riemann_HLLC(ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i], γ)
            elseif user_input.Solver_Input.riemanntype == :Exact
                f1, f2, f3 = ExactRiemannSolve!((ρR[i-1], uR[i-1], pR[i-1], ρL[i], uL[i], pL[i]),γ)
            end
            F1[i] = f1
            F2[i] = f2
            F3[i] = f3

        end
    end

    # Update conserved variables (interior cells)
    Threads.@threads for i in 2:nx-1
        @inbounds begin
            U.density_centers[i] -= dt/dx * (F1[i+1] - F1[i])
            U.momentum_centers[i] -= dt/dx * (F2[i+1] - F2[i])
            U.total_energy_centers[i] -= dt/dx * (F3[i+1] - F3[i])

        end
    end
end