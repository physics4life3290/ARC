




function LaxFriedrichs_Step(U::ConservativeVars, F::FluxVars, dt::Float64, ghost_zones::Int, total_zones::Int, spacing::Float64)
    for i in ghost_zones+1:total_zones-ghost_zones
        U.density[i] = 0.5 * (U.density[i+1] + U.density[i-1]) + dt/(2*spacing) * (F.dens_flux[i+1] - F.dens_flux[i-1])
        U.momentum[i] = 0.5 * (U.momentum[i+1] + U.momentum[i-1]) + dt/(2*spacing) * (F.mome_flux[i+1] - F.mome_flux[i-1])
        U.total_energy[i] = 0.5 * (U.total_energy[i+1] + U.total_energy[i-1]) + dt/(2*spacing) * (F.tot_ener_flux[i+1] - F.tot_ener_flux[i-1])
    end
end