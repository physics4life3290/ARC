




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


# New way of doing this via FLASH Methods paper.
function steepening_contact!(ρ::Vector{Float64}, P::Vector{Float64},
                             coord::Vector{Float64}, γ::Float64,
                             varL::Vector{Float64}, varR::Vector{Float64})

    N = length(ρ)

    # -----------------------------
    # 1. Compute face coords and Δξ
    # -----------------------------
    ξf = zeros(N+1)
    for i in 1:N-1
        ξf[i+1] = 0.5*(coord[i] + coord[i+1])
    end
    ξf[1]   = coord[1]   - 0.5*(coord[2]-coord[1])
    ξf[end] = coord[end] + 0.5*(coord[end]-coord[end-1])

    Δξ = similar(ρ)
    Δξ = [ξf[i+1] - ξf[i] for i in 1:N]

    # -----------------------------
    # 2. Loop over interior cells
    # -----------------------------
    Threads.@threads for i in 3:N-3   # need neighbors up to i±3
        @inbounds begin
            # Monotonized slopes
            δvar_im1 = slope_limiter(ρ[i-1] - ρ[i-2], ρ[i-1], ρ[i-2], ρ[i])
            δvar_ip1 = slope_limiter(ρ[i+1] - ρ[i], ρ[i+1], ρ[i], ρ[i+2])

            # Linear steepened interfaces
            ρL_d = ρ[i-1] + 0.5*δvar_im1
            ρR_d = ρ[i+1] - 0.5*δvar_ip1

            ξcumsum = cumsum(Δξ)  # once
            Δξ_im1 = ξcumsum[i-1]
            Δξ_i   = ξcumsum[i]
            Δξ_ip1 = ξcumsum[i+1]


            δ2ρ_i = (1 / (Δξ_im1 + Δξ_i + Δξ_ip1)) *
                    ((ρ[i+1] - ρ[i]) / (Δξ_ip1 + Δξ_i) -
                    (ρ[i] - ρ[i-1]) / (Δξ_i + Δξ_im1))

            # Neighbor second derivatives
            δ2ρ_im1 = (1 / (Δξ[i-2] + Δξ[i-1] + Δξ[i])) *
                    ((ρ[i] - ρ[i-1]) / (Δξ[i] + Δξ[i-1]) -
                    (ρ[i-1] - ρ[i-2]) / (Δξ[i-1] + Δξ[i-2]))

            δ2ρ_ip1 = (1 / (Δξ_ip1 + Δξ[i+2] + Δξ[i+3])) *
                    ((ρ[i+2] - ρ[i+1]) / (Δξ[i+3] + Δξ[i+2]) -
                    (ρ[i+1] - ρ[i])   / (Δξ[i+2] + Δξ_ip1))

            # Cell center positions
            ξ_im1 = sum(Δξ[1:i-1])
            ξ_i   = sum(Δξ[1:i])
            ξ_ip1 = sum(Δξ[1:i+1])

            # η̅_i
            ηbar = ((δ2ρ_ip1 - δ2ρ_im1) / (ξ_ip1 - ξ_im1)) *
                (((ξ_i - ξ_im1)^3 + (ξ_ip1 - ξ_i)^3) / (ρ[i+1] - ρ[i-1]))

            # Conditions
            if abs(ρ[i+1]-ρ[i-1]) / min(ρ[i+1],ρ[i-1]) < 0.01
                ηbar = 0.0
            elseif δ2ρ_im1*δ2ρ_ip1 > 0
                ηbar = 0.0
            elseif abs(P[i+1]-P[i-1]) / min(P[i+1],P[i-1]) >
                0.1*γ*abs(ρ[i+1]-ρ[i-1]) / min(ρ[i+1],ρ[i-1])
                ηbar = 0.0
            end

            # Steepening coefficient
            η_i = max(0.0, min(20*(ηbar - 0.05), 1.0))

            # Update interfaces
            varL[i] = varL[i]*(1 - η_i) + ρL_d*η_i
            varR[i] = varR[i]*(1 - η_i) + ρR_d*η_i
        end
    end

    return varL, varR
end
