




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