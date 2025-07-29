




# Solve star region
function solve_star_region(rhoL, uL, pL, rhoR, uR, pR, γ; max_iter=1E4)
    p_old = guess_pressure(rhoL, uL, pL, rhoR, uR, pR, γ)
    for iter in 1:max_iter
        fL, dFL = pressure_function(p_old, rhoL, pL, γ)
        fR, dFR = pressure_function(p_old, rhoR, pR, γ)
        f = fL + fR + (uR - uL)
        df = dFL + dFR
        p_new = p_old - f/df
        if abs(p_new - p_old) < 1e-10
            p_old = p_new
            break
        end
        p_old = max(1e-10, p_new)
    end
    p_star = p_old
    fL, _ = pressure_function(p_star, rhoL, pL, γ)
    fR, _ = pressure_function(p_star, rhoR, pR, γ)
    u_star = 0.5*(uL + uR + fR - fL)
    return p_star, u_star
end