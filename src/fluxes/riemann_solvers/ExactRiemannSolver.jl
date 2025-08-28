

#using Plots
include("../../solvers/Iterative/NewtonRaphson.jl")

cs(γ, p, ρ) = sqrt(γ * p / ρ)

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

# ---------------------------------------
# Pressure Estimates for Riemann Solvers
# ---------------------------------------

# Initial guess (PVRS)
function PVRS_guess(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    aL, aR = cs(γ, pL, ρL), cs(γ, pR, ρR)
    p̃ = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR) * (aL + aR)
    p0 = max(1e-8, p̃)

    return p0

end

# Two-Rarefaction Approximate Pressure (Toro Eq. 4.44)
function two_rarefaction_pressure(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    numerator = cL + cR - 0.5 * (γ - 1) * (uR - uL)
    denominator = cL / pL^((γ - 1) / (2γ)) + cR / pR^((γ - 1) / (2γ))

    p_pv = (numerator / denominator)^(2γ / (γ - 1))

    return max(1e-8, p_pv) # positivity safeguard

end

# Two-Shock Approximate Pressure (Toro Eq. 4.47)
function two_shock_pressure(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64)

    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    # Shock impedance factors (Eq. 4.47)
    gL = sqrt((2 / ((γ + 1) * ρL)) / (pL + (γ - 1) / (γ + 1) * pL))
    gR = sqrt((2 / ((γ + 1) * ρR)) / (pR + (γ - 1) / (γ + 1) * pR))

    numerator = gL * pL + gR * pR - (uR - uL)
    denominator = gL + gR

    p_ts = numerator / denominator

    return max(1e-8, p_ts) # positivity safeguard
    
end

function solve_star_pressure(ρL::Float64, uL::Float64, pL::Float64, ρR::Float64, uR::Float64, pR::Float64, γ::Float64; tol=1e-12, maxiter=10000)

    p0 = PVRS_guess(ρL, uL, pL, ρR, uR, pR, γ)

    f(p) = (wave_function(p, (ρL, uL, pL), γ)[1] + wave_function(p, (ρR, uR, pR), γ)[1] + (uR - uL))
    df(p) = (wave_function(p, (ρL, uL, pL), γ)[2] + wave_function(p, (ρR, uR, pR), γ)[2])

    pstar, converged, niter = newton_raphson(f, df, p0; tol=tol, maxiter=maxiter)
    
    if !converged
        @warn "Negative star pressure p* = $pstar, trying another guess"
        p0 = two_rarefaction_pressure(ρL, uL, pL, ρR, uR, pR, γ)
        pstar, converged, niter = newton_raphson(f, df, p0; tol=tol, maxiter=maxiter)
        if !converged
            @warn "Negative star pressure p* = $pstar, trying final guess"
            p0 = two_shock_pressure(ρL, uL, pL, ρR, uR, pR, γ)
            pstar, converged, niter = newton_raphson(f, df, p0; tol=tol, maxiter=maxiter)
            if !converged
                error("Newton–Raphson failed to converge after multiple attempts with p* = $pstar")
            end
        else
            error("Negative star pressure p* = $pstar, unable to find valid solution check input...")
        end
    end


    # Compute velocity in the star region
    fL, _ = wave_function(pstar, (ρL, uL, pL), γ)
    fR, _ = wave_function(pstar, (ρR, uR, pR), γ)
    ustar = 0.5 * (uL + uR) + 0.5 * (fR - fL)

    return pstar, ustar

end

function star_density(state::Tuple{Float64, Float64, Float64}, pstar::Float64, γ::Float64)
    ρ, u, p = state

    if pstar > p
        # Shock
        num = pstar / p + (γ - 1.0) / (γ + 1.0)
        den = ( (γ - 1.0) / (γ + 1.0) ) * (pstar / p) + 1.0
        return ρ * (num / den)
    else
        # Rarefaction
        return ρ * (pstar / p)^(1.0 / γ)
    end

end

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