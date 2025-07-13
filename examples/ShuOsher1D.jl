




function Construct1DShuOsherPrimitives(x::AbstractVector{<:Real}, γ::Real)
    total_zones = length(x)

    W = PrimitiveVariables(
        zeros(total_zones),   # density
        zeros(total_zones),   # velocity
        zeros(total_zones),   # pressure
        zeros(total_zones)    # internal energy
    )

    for i in 1:total_zones
        if x[i] < -4.0
            W.density[i] = 3.857143
            W.velocity[i] = 2.629369
            W.pressure[i] = 10.3333
        else
            W.density[i] = 1.0 + 0.2 * sin(5 * x[i])
            W.velocity[i] = 0.0
            W.pressure[i] = 1.0
        end
        W.int_energy[i] = W.pressure[i] / ((γ - 1) * W.density[i])
    end

    return W
end