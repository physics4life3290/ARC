# ---------------------------
# Utility Functions
# ---------------------------

function newton_raphson(f, df, x0; tol=1e-8, maxiter=1000)
    x = x0
    for _ in 1:maxiter
        fx, dfx = f(x), df(x)
        abs(dfx) < eps() && error("Derivative too close to zero.")
        dx = fx / dfx
        x -= dx
        abs(dx) < tol && return x
    end
    error("Newton-Raphson did not converge.")
end

function sound_speed(γ, p, ρ)
    sqrt(γ * p / ρ)
end

function shock_coefficients(ρ, p, γ)
    A = 2 / ((γ+1) * ρ)
    B = (γ-1) / (γ+1) * p
    return A, B
end

# ---------------------------
# Pressure Guess and Functions
# ---------------------------

function pressure_function(p_star, states, γ)
    ρL, uL, pL, ρR, uR, pR = states

    f(p, ρ, u, p₀) = begin
        if p > p₀
            A, B = shock_coefficients(ρ, p₀, γ)
            return (p - p₀) * sqrt(A / (p + B))
        else
            return (2 * sound_speed(γ, p₀, ρ) / (γ - 1)) * ((p / p₀)^((γ - 1) / (2 * γ)) - 1)
        end
    end

    return f(p_star, ρL, uL, pL) + f(p_star, ρR, uR, pR) + uR - uL
end


function pressure_derivative(p_star, states, γ)
    ρL, uL, pL, ρR, uR, pR = states

    df(p, ρ, p₀) = p > p₀ ? begin
            A, B = shock_coefficients(ρ, p₀, γ)
            (1 - (p - p₀) / (2 * (B + p))) * sqrt(A / (p + B))
        end :
        (1 / (ρ * sound_speed(γ, p₀, ρ))) * (p / p₀) ^ (-(γ + 1) / (2 * γ))

    return df(p_star, ρL, pL) + df(p_star, ρR, pR)
end


# ---------------------------
# Main Solver
# ---------------------------

function ExactRiemannSolve!(states, γ)

    ρL, uL, pL, ρR, uR, pR = states
    p_star_guess = 0.5 * (pL + pR)

    p_star = newton_raphson(
        p -> pressure_function(p, states, γ),
        p -> pressure_derivative(p, states, γ),
        p_star_guess
    )

    # For left side
    A_L, B_L = shock_coefficients(ρL, pL, γ)
    fL = p_star > pL ? 
        (p_star - pL) * sqrt(A_L / (p_star + B_L)) :
        (2 * sound_speed(γ, pL, ρL) / (γ - 1)) * ((p_star / pL)^((γ - 1) / (2 * γ)) - 1)

    # For right side
    A_R, B_R = shock_coefficients(ρR, pR, γ)
    fR = p_star > pR ? 
        (p_star - pR) * sqrt(A_R / (p_star + B_R)) :
        (2 * sound_speed(γ, pR, ρR) / (γ - 1)) * ((p_star / pR)^((γ - 1) / (2 * γ)) - 1)


    u_star = 0.5 * ((uL + uR) + (fR - fL))

    # Determine sampling region and return fluxes
    if u_star > 0.0
        if p_star > pL
            sL = uL - sound_speed(γ, pL, ρL) * sqrt((γ+1)/(2γ) * (p_star/pL) + (γ-1)/(2γ))
            if sL > 0
                ρ, u, p = ρL, uL, pL
            else
                ρ = ρL * ((p_star/pL + (γ-1)/(γ+1)) / ((γ-1)/(γ+1)*(p_star/pL) + 1))
                u, p = u_star, p_star
            end
        else
            cL = sound_speed(γ, pL, ρL)
            shL = uL - cL
            stL = u_star - cL * (p_star/pL)^((γ-1)/(2γ))
            if shL > 0
                ρ, u, p = ρL, uL, pL
            elseif stL < 0
                ρ = ρL * (p_star/pL)^(1/γ)
                u, p = u_star, p_star
            else
                ρ = ρL * ((2/(γ+1)) + (γ-1)/((γ+1)*cL) * (uL - u_star))^(2/(γ-1))
                u = (2/(γ+1)) * (cL + (γ-1)/2 * uL)
                p = pL * ((2/(γ+1)) + (γ-1)/((γ+1)*cL) * (uL - u_star))^(2γ/(γ-1))
            end
        end
    else
        if p_star > pR
            sR = uR + sound_speed(γ, pR, ρR) * sqrt((γ+1)/(2γ) * (p_star/pR) + (γ-1)/(2γ))
            if sR < 0
                ρ, u, p = ρR, uR, pR
            else
                ρ = ρR * ((p_star/pR + (γ-1)/(γ+1)) / ((γ-1)/(γ+1)*(p_star/pR) + 1))
                u, p = u_star, p_star
            end
        else
            cR = sound_speed(γ, pR, ρR)
            shR = uR + cR
            stR = u_star + cR * (p_star/pR)^((γ-1)/(2γ))
            if shR < 0
                ρ, u, p = ρR, uR, pR
            elseif stR > 0
                ρ = ρR * (p_star/pR)^(1/γ)
                u, p = u_star, p_star
            else
                ρ = ρR * ((2/(γ+1)) - (γ-1)/((γ+1)*cR) * (uR - u_star))^(2/(γ-1))
                u = (2/(γ+1)) * (-cR + (γ-1)/2 * uR)
                p = pR * ((2/(γ+1)) - (γ-1)/((γ+1)*cR) * (uR - u_star))^(2γ/(γ-1))
            end
        end
    end

    # Compute fluxes
    ρu = ρ * u
    momentum_flux = ρ * u^2 + p
    E = p / (γ - 1) + 0.5 * ρ * u^2
    energy_flux = u * (E + p)

    return ρu, momentum_flux, energy_flux
end
