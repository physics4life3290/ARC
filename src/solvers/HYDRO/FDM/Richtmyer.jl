




function RichtmyerStep!(W, U::ConservativeVariables,
                        F::FluxVariables,
                        _grid::CartesianGrid1D,
                        UserInput,
                        dt::Float64,
                        ghost_zones, total_zones, spacing, γ)

    CalculateFlux!(W, U,F)
    
    # Allocate U_half as a separate copy to avoid overwriting U
    U_half = ConservativeVariables(
        copy(U.density_centers),
        copy(U.momentum_centers),
        copy(U.total_energy_centers),
        nothing,
        nothing, 
        nothing
    )

    # Predictor step
    Threads.@threads for i in (ghost_zones + 1):(total_zones - ghost_zones)
        @inbounds begin
            ic, il = i, i-1 

            U_half.density_centers[i] = 0.5*(U.density_centers[ic] + U.density_centers[il]) -
                                        (dt/(2*spacing))*(F.density_flux[ic] - F.density_flux[il])
            U_half.momentum_centers[i] = 0.5*(U.momentum_centers[ic] + U.momentum_centers[il]) -
                                        (dt/(2*spacing))*(F.momentum_flux[ic] - F.momentum_flux[il])
            U_half.total_energy_centers[i] = 0.5*(U.total_energy_centers[ic] + U.total_energy_centers[il]) -
                                            (dt/(2*spacing))*(F.total_energy_flux[ic] - F.total_energy_flux[il])
        end
    end

    # Reconstruct primitive variables from U_half
    W_half = PrimitiveVariables(
        copy(U_half.density_centers),
        U_half.momentum_centers ./ U_half.density_centers,
        zeros(length(U_half.density_centers)),
        zeros(length(U_half.density_centers)),
        nothing, nothing, nothing, nothing
    )
    W_half.pressure_centers .= (γ - 1) .* (U_half.total_energy_centers .- 0.5 .* U_half.density_centers .* W_half.velocity_centers .^ 2)
    W_half.internal_energy_centers .= W_half.pressure_centers ./ ((γ-1) .* U_half.density_centers)

    apply_boundary_conditions(UserInput, U_half, _grid)

    # Compute F_half from U_half
    F_half = FluxVariables(
        U_half.momentum_centers,
        U_half.momentum_centers .* W_half.velocity_centers .+ W_half.pressure_centers,
        W_half.velocity_centers .* (U_half.total_energy_centers .+ W_half.pressure_centers)
    )
   

    # Corrector step
    Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
        @inbounds begin
            U.density_centers[i]      -= (dt/spacing)*(F_half.density_flux[i+1] - F_half.density_flux[i])
            U.momentum_centers[i]     -= (dt/spacing)*(F_half.momentum_flux[i+1] - F_half.momentum_flux[i])
            U.total_energy_centers[i] -= (dt/spacing)*(F_half.total_energy_flux[i+1] - F_half.total_energy_flux[i])
        end
    end

    # Final boundary conditions on updated U
    apply_boundary_conditions(UserInput, U, _grid)
end

