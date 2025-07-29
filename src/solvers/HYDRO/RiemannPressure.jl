




# Function to compute the pressure estimate using the two-rarefaction or two-shock approximation (PVRS)
function guess_pressure(rhoL, uL, pL, rhoR, uR, pR, γ)
    cL = sqrt(γ * pL / rhoL)
    cR = sqrt(γ * pR / rhoR)
    p_pv = 0.5*(pL + pR) - 0.125*(uR - uL)*(rhoL + rhoR)*(cL + cR)
    return max(1e-6, p_pv)
end

# Function f(p) and its derivative for Newton iteration
function pressure_function(p, rho, p_k, γ)
    if p > p_k
        A = 2/((γ+1)*rho)
        B = (γ-1)/(γ+1)*p_k
        f = (p - p_k)*sqrt(A/(p + B))
        df = sqrt(A/(p + B))*(1 - 0.5*(p - p_k)/(p + B))
    else
        a_k = sqrt(γ*p_k/rho)
        f = (2*a_k/(γ-1))*((p/p_k)^((γ-1)/(2γ)) - 1)
        df = (1/(rho*a_k))*(p/p_k)^(-(γ+1)/(2γ))
    end
    return f, df
end
