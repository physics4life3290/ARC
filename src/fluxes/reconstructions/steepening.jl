




function steepening(ρ::AbstractVector, p::AbstractVector, slope::AbstractVector)
    n = length(ρ)
    Threads.@threads for i in 2:(n-1)
        @inbounds begin
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
    end
    
    return slope
end

function compute_steepening_coefficient(density::AbstractVector, pressure::AbstractVector; ε=1e-10)
    n = length(density)
    β = zeros(n)  # start with no steepening

    Threads.@threads for i in 2:(n-1)
        @inbounds begin
            Δρ = abs(density[i+1] - density[i-1])
            Δp = abs(pressure[i+1] - pressure[i-1])
            ρ_avg = (density[i+1] + 2*density[i] + density[i-1]) / 4
            p_avg = (pressure[i+1] + 2*pressure[i] + pressure[i-1]) / 4

            # Normalize gradients
            grad_ρ = Δρ / (ρ_avg + ε)
            grad_p = Δp / (p_avg + ε)

            # Heuristic for contact: large density gradient but small pressure gradient
            if grad_ρ > 0.1 && grad_p < 0.05
                β[i] = 1.0
            elseif grad_ρ > 0.05 && grad_p < 0.1
                β[i] = 0.5
            else
                β[i] = 0.0
            end
        end
    end

    # Boundaries: no steepening
    β[1] = β[2]
    β[end] = β[end-1]

    return β
end
