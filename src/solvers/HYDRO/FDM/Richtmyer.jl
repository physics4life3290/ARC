




function RichtmyerStep!(U::ConservativeVars,
                        F::FluxVars,
                        dt::Float64,
                        ghost_zones, total_zones, spacing, γ)

    # allocate half–step (face-centered) flux arrays
    Nface = total_zones + 1
    dens_flux_half  = zeros(Nface)
    mome_flux_half  = zeros(Nface)
    ener_flux_half  = zeros(Nface)

    # Predictor: loop over faces i = ghost_zones … total_zones-ghost_zones
    for i in ghost_zones:(total_zones - ghost_zones + 1)
        # convert face index to cell center offsets
        ic   = i      # right cell
        il   = i - 1  # left  cell

        # compute U_half at the face
        ρ_half = 0.5*(U.density[ic] + U.density[il]) -
                 (dt/(2*spacing))*(F.dens_flux[ic] - F.dens_flux[il])
        m_half = 0.5*(U.momentum[ic] + U.momentum[il]) -
                 (dt/(2*spacing))*(F.mome_flux[ic] - F.mome_flux[il])
        E_half = 0.5*(U.total_energy[ic] + U.total_energy[il]) -
                 (dt/(2*spacing))*(F.tot_ener_flux[ic] - F.tot_ener_flux[il])

        u_half = m_half / ρ_half
        p_half = (γ - 1)*(E_half - 0.5*ρ_half*u_half^2)

        # store half-step flux at face i
        dens_flux_half[i] = m_half
        mome_flux_half[i] = m_half*u_half + p_half
        ener_flux_half[i] = u_half*(E_half + p_half)
    end

    # Corrector: update cells i = ghost_zones+1 … total_zones-ghost_zones
    for i in (ghost_zones+1):(total_zones-ghost_zones)
        U.density[i]      -= (dt/spacing)*(dens_flux_half[i+1] - dens_flux_half[i])
        U.momentum[i]     -= (dt/spacing)*(mome_flux_half[i+1] - mome_flux_half[i])
        U.total_energy[i] -= (dt/spacing)*(ener_flux_half[i+1] - ener_flux_half[i])
    end

end

