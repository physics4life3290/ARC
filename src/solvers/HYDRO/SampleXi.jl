




# Sample solution at xi = x/t
function sample_xi(xi, rhoL, uL, pL, rhoR, uR, pR, p_star, u_star, γ)
    if xi <= u_star
        if p_star > pL
            sL = uL - sqrt((γ+1)/(2γ)*p_star/pL + (γ-1)/(2γ))*sqrt(γ*pL/rhoL)
            return xi <= sL ? (rhoL, uL, pL) : (rhoL*((p_star/pL + (γ-1)/(γ+1))/((γ-1)/(γ+1)*p_star/pL + 1)), u_star, p_star)
        else
            aL = sqrt(γ*pL/rhoL)
            shL = uL - aL
            stL = u_star - sqrt(γ*p_star/rhoL)
            if xi <= shL
                return rhoL, uL, pL
            elseif xi <= stL
                u = (2/(γ+1))*(aL + (γ-1)/2*uL + xi)
                a = (2/(γ+1))*(aL + (γ-1)/2*(uL - xi))
                rho = rhoL*(a/aL)^(2/(γ-1))
                p = pL*(a/aL)^(2*γ/(γ-1))
                return rho, u, p
            else
                return rhoL*(p_star/pL)^(1/γ), u_star, p_star
            end
        end
    else
        if p_star > pR
            sR = uR + sqrt((γ+1)/(2γ)*p_star/pR + (γ-1)/(2γ))*sqrt(γ*pR/rhoR)
            return xi >= sR ? (rhoR, uR, pR) : (rhoR*((p_star/pR + (γ-1)/(γ+1))/((γ-1)/(γ+1)*p_star/pR + 1)), u_star, p_star)
        else
            aR = sqrt(γ*pR/rhoR)
            shR = uR + aR
            stR = u_star + sqrt(γ*p_star/rhoR)
            if xi >= shR
                return rhoR, uR, pR
            elseif xi >= stR
                u = (2/(γ+1))*(-aR + (γ-1)/2*uR + xi)
                a = (2/(γ+1))*(aR - (γ-1)/2*(uR - xi))
                rho = rhoR*(a/aR)^(2/(γ-1))
                p = pR*(a/aR)^(2*γ/(γ-1))
                return rho, u, p
            else
                return rhoR*(p_star/pR)^(1/γ), u_star, p_star
            end
        end
    end
end
