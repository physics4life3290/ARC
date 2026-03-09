# ============================================================================= #
#                                                                               #
#   HLLC Riemann Solver — Algorithms for Realistic Computations (ARC)          #
#   Toro (2009) Chapter 10 formulation                                          #
#                                                                               #
# ============================================================================= #

# ============================================================================= #
#   PRIMITIVE / CONSERVED / FLUX LAYER                                          #
#   Dispatches on AbstractEOS — defined in EOS module, not here                #
# ============================================================================= #

@inline function conserved(eos::AbstractEOS, ρ::Float64, u::Float64, p::Float64)
    E = total_energy(eos, ρ, u, p)
    return (ρ, ρ * u, E)
end

@inline function flux(eos::AbstractEOS, ρ::Float64, u::Float64, p::Float64)
    E = total_energy(eos, ρ, u, p)
    return (ρ * u, ρ * u^2 + p, u * (E + p))
end


# ============================================================================= #
#   WAVE SPEED ESTIMATES                                                         #
# ============================================================================= #

# Davis (1988) — simple, cheap, reliable for most flows
@inline function wavespeeds_Davis(eos::AbstractEOS,
                                   ρL::Float64, uL::Float64, pL::Float64,
                                   ρR::Float64, uR::Float64, pR::Float64)
    eL = specific_energy(eos, ρL, pL)
    eR = specific_energy(eos, ρR, pR)
    cL = sound_speed(eos, ρL, eL)
    cR = sound_speed(eos, ρR, eR)
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)
    return SL, SR
end

# Einfeldt (1988) — Roe-averaged, more robust near vacuum and strong rarefactions
@inline function wavespeeds_Einfeldt(eos::AbstractEOS,
                                      ρL::Float64, uL::Float64, pL::Float64,
                                      ρR::Float64, uR::Float64, pR::Float64)
    eL = specific_energy(eos, ρL, pL)
    eR = specific_energy(eos, ρR, pR)
    cL = sound_speed(eos, ρL, eL)
    cR = sound_speed(eos, ρR, eR)

    # Roe averages
    sqρL  = sqrt(ρL);  sqρR = sqrt(ρR)
    u_roe = (sqρL * uL + sqρR * uR) / (sqρL + sqρR)
    HL    = (total_energy(eos, ρL, uL, pL) + pL) / ρL
    HR    = (total_energy(eos, ρR, uR, pR) + pR) / ρR
    H_roe = (sqρL * HL + sqρR * HR) / (sqρL + sqρR)

    # Roe-averaged sound speed — generalized via enthalpy
    c_roe = sqrt(max(0.0, H_roe - 0.5 * u_roe^2))

    SL = min(uL - cL, u_roe - c_roe)
    SR = max(uR + cR, u_roe + c_roe)
    return SL, SR
end


# ============================================================================= #
#   STAR REGION                                                                  #
# ============================================================================= #

# Contact wave speed S* — Toro Eq. 10.58
@inline function S_star(ρL::Float64, uL::Float64, pL::Float64, SL::Float64,
                         ρR::Float64, uR::Float64, pR::Float64, SR::Float64)
    num = pR - pL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)
    den = ρL * (SL - uL) - ρR * (SR - uR)
    return num / den
end

# Star region conserved state — Toro Eq. 10.53
@inline function U_star_state(eos::AbstractEOS,
                               ρ::Float64, u::Float64, p::Float64,
                               S::Float64, Ss::Float64)
    ρ_star   = ρ * (S - u) / (S - Ss)
    mom_star = ρ_star * Ss
    E        = total_energy(eos, ρ, u, p)
    E_star   = ρ_star * (E / ρ + (Ss - u) * (Ss + p / (ρ * (S - u))))
    return (ρ_star, mom_star, E_star)
end

# Star region flux — Toro Eq. 10.71
@inline function F_star_flux(U::NTuple{3,Float64}, F::NTuple{3,Float64},
                              S::Float64, U_star::NTuple{3,Float64})
    return (
        F[1] + S * (U_star[1] - U[1]),
        F[2] + S * (U_star[2] - U[2]),
        F[3] + S * (U_star[3] - U[3])
    )
end


# ============================================================================= #
#   VACUUM GUARD                                                                 #
# ============================================================================= #

const VACUUM_THRESHOLD = 1.0e-14

@inline function check_vacuum(ρL::Float64, ρR::Float64)
    if ρL < VACUUM_THRESHOLD || ρR < VACUUM_THRESHOLD
        error("HLLC: vacuum state detected (ρL=$ρL, ρR=$ρR). Vacuum Riemann solver not yet implemented.")
    end
end


# ============================================================================= #
#   HLLC RIEMANN SOLVER                                                          #
#   Toro (2009) Chapter 10                                                       #
#                                                                               #
#   Arguments:                                                                  #
#     eos          — AbstractEOS instance (IdealGasEOS, StiffenedGasEOS, ...)   #
#     ρL, uL, pL  — left primitive state                                        #
#     ρR, uR, pR  — right primitive state                                       #
#     wavespeed    — :Davis (default) or :Einfeldt                              #
#                                                                               #
#   Returns:                                                                    #
#     HLLC numerical flux (ρ, ρu, E) at the interface                          #
# ============================================================================= #

function Riemann_HLLC(eos::AbstractEOS,
                      ρL::Float64, uL::Float64, pL::Float64,
                      ρR::Float64, uR::Float64, pR::Float64;
                      wavespeed::Symbol = :Einfeldt)

    check_vacuum(ρL, ρR)

    # --- Wave speed estimates ---
    SL, SR = if wavespeed === :Einfeldt
        wavespeeds_Einfeldt(eos, ρL, uL, pL, ρR, uR, pR)
    else
        wavespeeds_Davis(eos, ρL, uL, pL, ρR, uR, pR)
    end

    # --- Contact wave speed ---
    Ss = S_star(ρL, uL, pL, SL, ρR, uR, pR, SR)

    # --- Conserved states and fluxes ---
    UL = conserved(eos, ρL, uL, pL)
    UR = conserved(eos, ρR, uR, pR)
    FL = flux(eos, ρL, uL, pL)
    FR = flux(eos, ρR, uR, pR)

    # --- HLLC flux selection ---
    if SL >= 0.0
        return FL

    elseif SL <= 0.0 <= Ss
        UL_star = U_star_state(eos, ρL, uL, pL, SL, Ss)
        return F_star_flux(UL, FL, SL, UL_star)

    elseif Ss <= 0.0 <= SR
        UR_star = U_star_state(eos, ρR, uR, pR, SR, Ss)
        return F_star_flux(UR, FR, SR, UR_star)

    else
        return FR
    end
end