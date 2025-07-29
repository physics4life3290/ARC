




function LaxFriedrichs_Step(U::ConservativeVariables, F::FluxVariables, dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64)

    Threads.@threads for i in ghost_zones+1:total_zones-ghost_zones
        @inbounds begin
            U.density_centers[i]      = 0.5*(U.density_centers[i+1] + U.density_centers[i-1]) -
                                (dt/(2*spacing))*(F.density_flux[i+1]   - F.density_flux[i-1])
            U.momentum_centers[i]     = 0.5*(U.momentum_centers[i+1] + U.momentum_centers[i-1]) -
                                (dt/(2*spacing))*(F.momentum_flux[i+1]   - F.momentum_flux[i-1])
            U.total_energy_centers[i] = 0.5*(U.total_energy_centers[i+1] + U.total_energy_centers[i-1]) -
                                (dt/(2*spacing))*(F.total_energy_flux[i+1] - F.total_energy_flux[i-1])
        end
    end
end