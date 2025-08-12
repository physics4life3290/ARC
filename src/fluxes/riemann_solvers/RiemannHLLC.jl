




function Riemann_HLLC(ρL, uL, pL, ρR, uR, pR, γ)
    # Compute sound speeds
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    # Conservative variables (ρ, ρu, E)
    function conserved(ρ, u, p)
        E = p/(γ-1) + 0.5 * ρ * u^2
        return (ρ, ρ*u, E)
    end

    # Flux vector (ρu, ρu^2 + p, u(E+p))
    function flux(ρ, u, p)
        E = p/(γ-1) + 0.5 * ρ * u^2
        return (ρ*u, ρ*u^2 + p, u*(E + p))
    end

    UL = conserved(ρL, uL, pL)
    UR = conserved(ρR, uR, pR)
    FL = flux(ρL, uL, pL)
    FR = flux(ρR, uR, pR)

    # Estimate wave speeds (e.g. Davis' estimates or Einfeldt's)
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    # Compute pressure estimate in star region (Toro Eq. 10.58)
    numerator = pR - pL + ρL*uL*(SL - uL) - ρR*uR*(SR - uR)
    denominator = ρL*(SL - uL) - ρR*(SR - uR)
    S_star = numerator / denominator

    # Compute star region conservative states for left and right
    function U_star(ρ, u, p, S, S_star)
        ρ_star = ρ * (S - u) / (S - S_star)
        mom_star = ρ_star * S_star
        E = p/(γ-1) + 0.5 * ρ * u^2
        E_star = ((S - u)*E - p*u + p * S_star) / (S - S_star)
        return (ρ_star, mom_star, E_star)
    end

    # Flux in star region
    function F_star(U, F, S, S_star, U_star)
        return (
            F[1] + S * (U_star[1] - U[1]),
            F[2] + S * (U_star[2] - U[2]),
            F[3] + S * (U_star[3] - U[3])
        )
    end

    if SL >= 0
        return FL
    elseif SL <= 0 <= S_star
        U_L_star = U_star(ρL, uL, pL, SL, S_star)
        return F_star(UL, FL, SL, S_star, U_L_star)
    elseif S_star <= 0 <= SR
        U_R_star = U_star(ρR, uR, pR, SR, S_star)
        return F_star(UR, FR, SR, S_star, U_R_star)
    else
        return FR
    end
end