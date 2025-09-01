




function wave_function(p::Float64, state::Tuple{Float64, Float64, Float64}, γ::Float64)

    ρ, u, p_i = state
    a = cs(γ, p_i, ρ)
     
    if p > p_i
        # Shock
        A = 2.0 / ((γ + 1.0) * ρ)
        B = (γ - 1.0) / (γ + 1.0) * p_i
        sqrt_term = sqrt(A / (p + B))
        f  = (p - p_i) * sqrt_term
        df = sqrt_term * (1.0 - 0.5 * (p - p_i) / (p + B))
        return f, df
    else
        # Rarefaction
        pr = p / p_i
        expo = (γ - 1.0) / (2.0 * γ)
        f  = (2.0 * a / (γ - 1.0)) * (pr^expo - 1.0)
        df = (1.0 / (ρ * a)) * pr^(-(γ + 1.0) / (2.0 * γ))
        return f, df
    end

end