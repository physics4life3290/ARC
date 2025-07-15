




function RichtmyerStep(U::ConservativeVars, F::FluxVars, dt::Float64, ghost_zones, total_zones, spacing, γ)
    for i in ghost_zones+1:total_zones-ghost_zones
        # Predictor Step: compute midpoint conservative variables
        ρ_half = 0.5 * (U.density[i+1] + U.density[i-1]) + (dt / (2 * spacing)) * (F.dens_flux[i+1] - F.dens_flux[i-1])
        m_half = 0.5 * (U.momentum[i+1] + U.momentum[i-1]) + (dt / (2 * spacing)) * (F.mome_flux[i+1] - F.mome_flux[i-1])
        E_half = 0.5 * (U.total_energy[i+1] + U.total_energy[i-1]) + (dt / (2 * spacing)) * (F.tot_ener_flux[i+1] - F.tot_ener_flux[i-1])

        # Compute primitive variables from midpoint conservative vars
        u_half = m_half / ρ_half
        kinetic_half = 0.5 * ρ_half * u_half^2
        p_half = (γ - 1) * (E_half - kinetic_half)

        # Compute fluxes at the midpoint
        flux_ρ = m_half
        flux_m = m_half * u_half + p_half
        flux_E = u_half * (E_half + p_half)

        F.dens_flux[i] = flux_ρ
        F.mome_flux[i] = flux_m 
        F.tot_ener_flux[i] = flux_E
    end
    
    for i in ghost_zones+1:total_zones-ghost_zones
    # Corrector Step: update U using midpoint flux
        U.density[i]  -= (dt / spacing) * (F.dens_flux[i+1] - F.dens_flux[i-1])
        U.momentum[i] -= (dt / spacing) * (F.mome_flux[i+1] - F.mome_flux[i-1])
        U.total_energy[i] -= (dt / spacing) * (F.tot_ener_flux[i+1] - F.tot_ener_flux[i-1])
    end
end
