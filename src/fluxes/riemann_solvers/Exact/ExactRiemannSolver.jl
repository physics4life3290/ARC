

#using Plots
include("../../../solvers/Iterative/NewtonRaphson.jl")

cs(γ, p, ρ) = sqrt(γ * p / ρ)
left_shock_speed(ρL, uL, pL, pstar, γ) = uL - cs(γ, pL, ρL) * sqrt( ( (γ + 1.0) / (2.0 * γ) ) * (pstar / pL) + (γ - 1.0) / (2.0 * γ) )
right_shock_speed(ρR, uR, pR, pstar, γ) = uR + cs(γ, pR, ρR) * sqrt( ( (γ + 1.0) / (2.0 * γ) ) * (pstar / pR) + (γ - 1.0) / (2.0 * γ) )

function sample_exact(xi::Float64, ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, pstar::Float64, ustar::Float64, γ::Float64)
    # Precompute star states
    ρLstar = star_density((ρL, uL, pL), pstar, γ)
    ρRstar = star_density((ρR, uR, pR), pstar, γ)
    aLstar = cs(γ, pstar, ρLstar)
    aRstar = cs(γ, pstar, ρRstar)

    if xi <= ustar
        # LEFT side of contact
        if pstar > pL
            # Left shock
            SL = left_shock_speed(ρL, uL, pL, pstar, γ)
            if xi < SL
                return (ρL, uL, pL)
            else
                return (ρLstar, ustar, pstar)
            end
        else
            # Left rarefaction
            SHL = uL - aLstar          # head
            STL = ustar - aLstar     # tail
            if xi < SHL
                return (ρL, uL, pL)
            elseif xi > STL
                return (ρLstar, ustar, pstar)
            else
                # Inside left fan (use Riemann invariants)
                u = (2.0 / (γ + 1.0)) * (aLstar + 0.5*(γ - 1.0)*uL + xi)
                a = (2.0 / (γ + 1.0)) * (aLstar + 0.5*(γ - 1.0)*(uL - xi))
                ρ = ρL * (a / aLstar)^(2.0 / (γ - 1.0))
                p = γ^-1 * a^2 * ρ
                return (ρ, u, p)
            end
        end
    else
        # RIGHT side of contact
        if pstar > pR
            # Right shock
            SR = right_shock_speed(ρR, uR, pR, pstar, γ)
            if xi > SR
                return (ρR, uR, pR)
            else
                return (ρRstar, ustar, pstar)
            end
        else
            # Right rarefaction
            SHR = uR + aRstar          # head
            STR = ustar + aRstar     # tail
            if xi > SHR
                return (ρR, uR, pR)
            elseif xi < STR
                return (ρRstar, ustar, pstar)
            else
                # Inside right fan
                u = (2.0 / (γ + 1.0)) * (-aRstar + 0.5*(γ - 1.0)*uR + xi)
                a = (2.0 / (γ + 1.0)) * (aRstar - 0.5*(γ - 1.0)*(uR - xi))
                ρ = ρR * (a / aRstar)^(2.0 / (γ - 1.0))
                p = γ^-1 * a^2 * ρ
                return (ρ, u, p)
            end
        end
    end
end