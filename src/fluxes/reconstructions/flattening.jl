




function flattening(W, slope_ρ)
    nx = length(W.density_centers)
    for i in 2:nx-1
        shock_strength = abs(W.density_centers[i+1] - 2*W.density_centers[i] + W.density_centers[i-1])  # curvature of density
        flattening_coef = clamp(1.0 - 10.0 * shock_strength, 0.0, 1.0)
        slope_ρ[i] = slope_ρ[i] * flattening_coef
    end
end