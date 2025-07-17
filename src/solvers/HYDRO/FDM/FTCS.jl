




#FTCS Method 
function FTCS_Step!(U::ConservativeVars,
                    F::FluxVars,
                    dt::Float64,
                    ghost_zones::Int,
                    total_zones::Int,
                    spacing::Float64)

    # 1) copy old solution
    U_old = deepcopy(U)

    # 2) compute new solution from old
    for i in (ghost_zones+1):(total_zones-ghost_zones)
        U.density[i]      = U_old.density[i]      - (dt/(2*spacing))*(F.dens_flux[i+1]  - F.dens_flux[i-1])
        U.momentum[i]     = U_old.momentum[i]     - (dt/(2*spacing))*(F.mome_flux[i+1]  - F.mome_flux[i-1])
        U.total_energy[i] = U_old.total_energy[i] - (dt/(2*spacing))*(F.tot_ener_flux[i+1] - F.tot_ener_flux[i-1])
    end

    return nothing
end
