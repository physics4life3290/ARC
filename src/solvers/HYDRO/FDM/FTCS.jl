




#FTCS Method 
function FTCS_Step!(U::ConservativeVariables,
                    F::FluxVariables,
                    dt::Float64,
                    ghost_zones::Int,
                    total_zones::Int,
                    spacing::Float64)

    # 1) copy old solution
    U_old = deepcopy(U)

    # 2) compute new solution from old
    Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
        @inbounds begin
            U.density_centers[i]      = U_old.density_centers[i]      - (dt/(2*spacing))*(F.density_flux[i+1]  - F.density_flux[i-1])
            U.momentum_centers[i]     = U_old.momentum_centers[i]     - (dt/(2*spacing))*(F.momentum_flux[i+1]  - F.momentum_flux[i-1])
            U.total_energy_centers[i] = U_old.total_energy_centers[i] - (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
        end
    end

    return nothing
end
