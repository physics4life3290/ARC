




function cubic_interpolation(var::Vector{Float64}, Δξ::Vector{Float64}, i::Int)
    # shorthand
    Δξm = Δξ[i-1]   # Δξ_(i-1)
    Δξi = Δξ[i]     # Δξ_i
    Δξip1 = Δξ[i+1] # Δξ_(i+1)
    Δξip2 = Δξ[i+2] # Δξ_(i+2)

    vi   = var[i]
    vip1 = var[i+1]
    vim1 = var[i-1]

    # slope δvar_i
    δvar_i = (Δξi / (Δξm + Δξi + Δξip1)) * (
        ((2Δξm + Δξi) / (Δξip1 + Δξi)) * (vip1 - vi) +
        ((Δξi + 2Δξm) / (Δξm + Δξi)) * (vi - vim1)
    )

    # slope δvar_(i+1)
    vip2 = var[i+2]
    δvar_ip1 = (Δξip1 / (Δξi + Δξip1 + Δξip2)) * (
        ((2Δξi + Δξip1) / (Δξip2 + Δξip1)) * (vip2 - vip1) +
        ((Δξip1 + 2Δξi) / (Δξi + Δξip1)) * (vip1 - vi)
    )

    # Z terms
    Z1 = (Δξm + Δξi) / (2Δξi + Δξip1)
    Z2 = (Δξip2 + Δξip1) / (2Δξip1 + Δξi)

    ΔX = Δξm + Δξi + Δξip1 + Δξip2

    # cubic interpolation for var_R
    var_R = vi +
        (Δξi / (Δξi + Δξip1)) * (vip1 - vi) +
        (1/ΔX) * (
            (2Δξip1*Δξi)/(Δξip1+Δξi) * (Z1 - Z2) * (vip1 - vi)
            - Δξi * Z1 * δvar_ip1
            + Δξip1 * Z2 * δvar_i
        )

    return var_R, δvar_i, δvar_ip1
end


function slope_limiter(δvar, vi, vim1, vip1)
    if (vip1 - vi) * (vi - vim1) > 0
        return sign(δvar) * min(abs(δvar), 2abs(vi - vim1), 2abs(vip1 - vi))
    else
        return 0.0
    end
end

function parabolic_reconstruction(varL, varR, vi, ξ, ξL, ξR)
    α = (ξ - ξL) / (ξR - ξL)
    Δvar = varR - varL
    var6i = 6 * (vi - 0.5*(varL + varR))
    return varL + α*(Δvar + var6i*(1 - α))
end

function reconstruct_interfaces(var::Vector{Float64}, coord::Vector{Float64})
    N = length(var)

    # -----------------------------
    # 1. Compute face coordinates
    # -----------------------------
    ξf = zeros(N+1) # N+1 faces for N cells
    Threads.@threads for i in 1:N-1
        @inbounds begin
            ξf[i+1] = 0.5*(coord[i] + coord[i+1])
        end
    end
    # extrapolate outermost faces (assuming uniform extension at boundaries)
    ξf[1]   = coord[1]   - 0.5*(coord[2]-coord[1])
    ξf[end] = coord[end] + 0.5*(coord[end]-coord[end-1])

    # -----------------------------
    # 2. Cell widths Δξ
    # -----------------------------
    Δξ = [ξf[i+1] - ξf[i] for i in 1:N]

    # -----------------------------
    # 3. Allocate interface arrays
    # -----------------------------
    varL = zeros(N)
    varR = zeros(N)

    # -----------------------------
    # 4. Reconstruction loop
    # -----------------------------
    Threads.@threads for i in 2:N-2
        @inbounds begin 
            vR, δvar_i, δvar_ip1 = cubic_interpolation(var, Δξ, i)

            # slope limiter on δvar_i
            δvar_i_lim = slope_limiter(δvar_i, var[i], var[i-1], var[i+1])

            # assign right interface of zone i
            varR[i] = vR

            # assign left interface of zone i+1 (mirror)
            varL[i+1] = varR[i]
        end
    end

    return varL, varR, Δξ, ξf
end


#using Plots

#= -----------------------------
# Define grid and function
# -----------------------------
N = 50
ξ = range(0, 2π, length=N)
Δξ = fill(ξ[2]-ξ[1], N)      # uniform spacing
var = sin.(3*ξ) + 0.3*sin.(10*ξ) # exotic oscillatory function

# -----------------------------
# Compute reconstructed interfaces
# -----------------------------
varR = zeros(N)
varL = zeros(N)

for i in 2:N-2
    # cubic interpolation for right interface
    vR, δvar_i, δvar_ip1 = cubic_interpolation(var, Δξ, i)
    # apply limiter
    δvar_i_lim = slope_limiter(δvar_i, var[i], var[i-1], var[i+1])
    varR[i] = vR

    # Left interface can be mirrored
    varL[i+1] = varR[i]
end

# -----------------------------
# Evaluate parabolic reconstruction
# -----------------------------
ξ_fine = range(0, 2π, length=500)
var_recon = zeros(length(ξ_fine))

for i in 2:N-2
    # cell boundaries
    ξL = ξ[i] - Δξ[i]/2
    ξR = ξ[i] + Δξ[i]/2

    # find fine grid points inside this cell
    idx = findall(x -> ξL <= x <= ξR, ξ_fine)
    for j in idx
        var_recon[j] = parabolic_reconstruction(varL[i], varR[i], var[i], ξ_fine[j], ξL, ξR)
    end
end

#-----------------------------
# Plot
# -----------------------------
plot(ξ, var, lw=2, marker=:circle, label="Original cell avg")
plot!(ξ_fine, var_recon, lw=2, label="Parabolic reconstruction", color=:red)
xlabel!("ξ")
ylabel!("var(ξ)")
title!("Exotic function cubic interpolation with parabolic reconstruction")
=#