




#FTCS Method 
function FTCS_Step(U::ConservativeVars, F::FluxVars, dt::Float64, ghost_zones, total_zones, spacing)
    for i in ghost_zones+1:total_zones-ghost_zones
        U.density[i] = U.density[i] + dt/(2*spacing) * (F.dens_flux[i+1] - F.dens_flux[i-1])
        U.momentum[i] = U.momentum[i] + dt/(2*spacing) * (F.mome_flux[i+1] - F.mome_flux[i-1])
        U.total_energy[i] = U.total_energy[i] + dt/(2*spacing) * (F.tot_ener_flux[i+1] - F.tot_ener_flux[i-1])
    end
end