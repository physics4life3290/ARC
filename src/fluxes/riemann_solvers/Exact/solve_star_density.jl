




function star_density(state::Tuple{Float64, Float64, Float64}, pstar::Float64, γ::Float64)
    ρ, u, p = state

    if pstar > p
        # Shock
        num = pstar / p + (γ - 1.0) / (γ + 1.0)
        den = ( (γ - 1.0) / (γ + 1.0) ) * (pstar / p) + 1.0
        return ρ * (num / den)
    else
        # Rarefaction
        return ρ * (pstar / p)^(1.0 / γ)
    end

end