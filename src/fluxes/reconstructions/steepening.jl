




function steepening(ρ, p, slope, iters)
    for i in 2:iters-1
        dCρ = 0.5 * (ρ[i+1] - ρ[i-1])
        ΔpL = abs(p[i] - p[i-1])
        ΔpR = abs(p[i+1] - p[i])
        ΔρL = abs(ρ[i]-ρ[i-1])
        ΔρR = abs(ρ[i+1]-ρ[i])
        steepening_coef = 0.0

        if (ΔρL > 0.02 * ρ[i]) && (ΔpL < 0.1 * p[i]) &&
            (ΔρR > 0.02 * ρ[i]) && (ΔpR < 0.1 * p[i])
            steepening_coef = 0.3  # strengthen slope slightly
        end
        slope[i] += steepening_coef * dCρ
    end
    
    return slope
end