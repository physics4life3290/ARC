




function Construct1DSedovBlastPrimitives(x::AbstractVector{<:Real}, dx::Real,
                                        blast_center::Real, blast_energy::Real,
                                        background_density::Real, background_pressure::Real,
                                        γ::Real; blast_width::Real = dx)

    total_zones = length(x)

    W = PrimitiveVariables(
        fill(background_density, total_zones),     # density
        zeros(total_zones),                        # velocity
        fill(background_pressure, total_zones),    # pressure
        zeros(total_zones)                         # internal energy
    )

    # Energy deposition over small region near blast center
    for i in 1:total_zones
        if abs(x[i] - blast_center) <= blast_width/2
            volume = dx
            W.pressure[i] = (γ - 1) * blast_energy / volume
        end
        W.int_energy[i] = W.pressure[i] / ((γ - 1) * W.density[i])
    end

    return W
end

