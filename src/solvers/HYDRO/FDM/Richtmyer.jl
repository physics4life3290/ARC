




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

function RichtmyerStep_AV!(W, U::ConservativeVariables,
                           F::FluxVariables,
                           _grid::CartesianGrid1D,
                           UserInput,
                           dt::Float64,
                           ghost_zones, total_zones, spacing, γ)

    # Compute fluxes at current time
    CalculateFlux!(W, U, F)

    # Characteristic "a" = max wavespeed (for viscosity coefficient)
    cs = sqrt.(γ .* W.pressure_centers ./ W.density_centers) # sound speed
    a  = maximum(abs.(W.velocity_centers) .+ cs)

    ν = 0.5 * (a * dt)^2 / spacing   # Richtmyer artificial viscosity coefficient

    Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
        @inbounds begin
            # Central flux difference (non-viscous term)
            U.density_centers[i]      -= (dt/(2*spacing)) * (F.density_flux[i+1] - F.density_flux[i-1])
            U.momentum_centers[i]     -= (dt/(2*spacing)) * (F.momentum_flux[i+1] - F.momentum_flux[i-1])
            U.total_energy_centers[i] -= (dt/(2*spacing)) * (F.total_energy_flux[i+1] - F.total_energy_flux[i-1])

            # Artificial viscosity term (second derivative form)
            U.density_centers[i]      += ν * (U.density_centers[i+1] - 2*U.density_centers[i] + U.density_centers[i-1]) / spacing^2
            U.momentum_centers[i]     += ν * (U.momentum_centers[i+1] - 2*U.momentum_centers[i] + U.momentum_centers[i-1]) / spacing^2
            U.total_energy_centers[i] += ν * (U.total_energy_centers[i+1] - 2*U.total_energy_centers[i] + U.total_energy_centers[i-1]) / spacing^2
        end
    end

    # Apply boundaries
    apply_boundary_conditions(UserInput, U, _grid)
end

function RichtmyerFVS!(W, U::ConservativeVariables,
                        _grid::CartesianGrid1D,
                        UserInput,
                        dt::Float64,
                        ghost_zones, total_zones, spacing, γ)

    # Compute F+ and F- from W
    F_plus, F_minus = FluxVectorSplit(W, γ)

    # Predictor step: compute U_half at i+1/2 using FVS upwinding
    U_half = ConservativeVariables(
        copy(U.density_centers),
        copy(U.momentum_centers),
        copy(U.total_energy_centers),
        nothing, nothing, nothing
    )

    Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
        @inbounds begin
            # Use FVS: right-moving flux uses i, left-moving flux uses i-1
            U_half.density_centers[i]      = 0.5*(U.density_centers[i] + U.density_centers[i-1]) -
                                              (dt/(2*spacing))*(F_plus.density_flux[i] - F_plus.density_flux[i-1] +
                                                                 F_minus.density_flux[i] - F_minus.density_flux[i-1])

            U_half.momentum_centers[i]     = 0.5*(U.momentum_centers[i] + U.momentum_centers[i-1]) -
                                              (dt/(2*spacing))*(F_plus.momentum_flux[i] - F_plus.momentum_flux[i-1] +
                                                                 F_minus.momentum_flux[i] - F_minus.momentum_flux[i-1])

            U_half.total_energy_centers[i] = 0.5*(U.total_energy_centers[i] + U.total_energy_centers[i-1]) -
                                              (dt/(2*spacing))*(F_plus.total_energy_flux[i] - F_plus.total_energy_flux[i-1] +
                                                                 F_minus.total_energy_flux[i] - F_minus.total_energy_flux[i-1])
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
    W_half.pressure_centers .= (γ - 1) .* (U_half.total_energy_centers .- 0.5 .* U_half.density_centers .* W_half.velocity_centers.^2)
    W_half.internal_energy_centers .= W_half.pressure_centers ./ ((γ-1) .* U_half.density_centers)

    apply_boundary_conditions(UserInput, U_half, _grid)

    # Compute F_half+ and F_half- from W_half
    Fp_half, Fm_half = FluxVectorSplit(W_half, γ)

    # Corrector step: update U using FVS upwinded differences
    Threads.@threads for i in (ghost_zones+1):(total_zones-ghost_zones)
        @inbounds begin
            U.density_centers[i]      -= (dt/spacing)*(Fp_half.density_flux[i+1] - Fp_half.density_flux[i] +
                                                       Fm_half.density_flux[i] - Fm_half.density_flux[i-1])
            U.momentum_centers[i]     -= (dt/spacing)*(Fp_half.momentum_flux[i+1] - Fp_half.momentum_flux[i] +
                                                       Fm_half.momentum_flux[i] - Fm_half.momentum_flux[i-1])
            U.total_energy_centers[i] -= (dt/spacing)*(Fp_half.total_energy_flux[i+1] - Fp_half.total_energy_flux[i] +
                                                       Fm_half.total_energy_flux[i] - Fm_half.total_energy_flux[i-1])
        end
    end

    apply_boundary_conditions(UserInput, U, _grid)
end
