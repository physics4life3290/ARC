




function Riemann_HLL(ρL::Float64, uL::Float64, pL::Float64,
                     ρR::Float64, uR::Float64, pR::Float64,
                     eos::AbstractEOS)

    # Toro's method - simplified pressure estimate
    p_star = max(1e-6, 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR))
    u_star = 0.5 * (uL + uR) + 0.5 * (pL - pR) / (ρL + ρR)

    # Wave speeds via EOS sound speed
    eL     = specific_energy(eos, ρL, pL)
    eR     = specific_energy(eos, ρR, pR)
    e_star = specific_energy(eos, ρL, p_star)  # uses ρL as representative density

    cL = sound_speed(eos, ρL, eL)
    cR = sound_speed(eos, ρR, eR)
    cL_star = sound_speed(eos, ρL, e_star)
    cR_star = sound_speed(eos, ρR, specific_energy(eos, ρR, p_star))

    SL = min(uL - cL, u_star - cL_star)
    SR = max(uR + cR, u_star + cR_star)

    if SL >= 0
        return flux(ρL, uL, pL, eos)
    elseif SR <= 0
        return flux(ρR, uR, pR, eos)
    else
        FL = flux(ρL, uL, pL, eos)
        FR = flux(ρR, uR, pR, eos)
        UL = conserved(ρL, uL, pL, eos)
        UR = conserved(ρR, uR, pR, eos)
        F1 = (SR * FL[1] - SL * FR[1] + SL * SR * (UR[1] - UL[1])) / (SR - SL)
        F2 = (SR * FL[2] - SL * FR[2] + SL * SR * (UR[2] - UL[2])) / (SR - SL)
        F3 = (SR * FL[3] - SL * FR[3] + SL * SR * (UR[3] - UL[3])) / (SR - SL)
        return F1, F2, F3
    end

end