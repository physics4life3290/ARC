





function Godunov_Step!(UserInput, _grid, W, U)
    γ = UserInput.secondary_input.γ
    dx = _grid.xcoord.spacing
    c = sqrt.(UserInput.secondary_input.γ .* W.pressure_centers ./ W.density_centers)
    dt = UserInput.secondary_input.cfl * _grid.xcoord.spacing / maximum(abs.(W.velocity_centers) .+ c)
    dt = min(dt, t_final - t)
    nx = _grid.xcoord.total_zones
    # Compute fluxes at interfaces
    F1 = zeros(nx+1)
    F2 = zeros(nx+1)
    F3 = zeros(nx+1)

    Threads.@threads for i in 2:nx
        @inbounds begin
            ρL, uL, pL = W.density_centers[i-1], W.velocity_centers[i-1], W.pressure_centers[i-1]
            ρR, uR, pR = W.density_centers[i], W.velocity_centers[i], W.pressure_centers[i]
            #f1, f2, f3 = Riemann_HLL(ρL, uL, pL, ρR, uR, pR)
            f1, f2, f3 = ExactRiemannSolve!((ρL, uL, pL, ρR, uR, pR), γ)
            #f1, f2, f3 = Riemann_HLLC(ρL, uL, pL, ρR, uR, pR, γ)
            F1[i] = f1
            F2[i] = f2
            F3[i] = f3
        end
    end

    # Update conserved variables
    Threads.@threads for i in 2:nx-1
        @inbounds begin
            U.density_centers[i] -= dt/dx * (F1[i+1] - F1[i])
            U.momentum_centers[i] -= dt/dx * (F2[i+1] - F2[i])
            U.total_energy_centers[i] -= dt/dx * (F3[i+1] - F3[i])
        end
    end
end