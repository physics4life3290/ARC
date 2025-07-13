




function Construct1DBlastWaveInteractionPrimitives(x::AbstractVector{<:Real}, γ::Real)
    total_zones = length(x)

    W = PrimitiveVariables(
        zeros(total_zones),   # density
        zeros(total_zones),   # velocity
        zeros(total_zones),   # pressure
        zeros(total_zones)    # internal energy
    )

    for i in 1:total_zones
        if x[i] < 0.1
            W.density[i] = 1.0
            W.velocity[i] = 0.0
            W.pressure[i] = 1000.0
        elseif x[i] < 0.9
            W.density[i] = 1.0
            W.velocity[i] = 0.0
            W.pressure[i] = 0.01
        else
            W.density[i] = 1.0
            W.velocity[i] = 0.0
            W.pressure[i] = 100.0
        end
        W.int_energy[i] = W.pressure[i] / ((γ - 1) * W.density[i])
    end

    return W
end