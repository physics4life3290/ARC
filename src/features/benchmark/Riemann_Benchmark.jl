




function Riemann_Benchmark(x::Vector, xc::Vector, ρ::Vector, u::Vector, p::Vector, γ::Float64, t::Float64)
    N = length(x)
    M = length(xc)


    ρ_sol = similar(x)
    u_sol = similar(x)
    p_sol = similar(x)

    # Loop over spatial points where we want solution
    @inbounds for j in 1:M
        xj = x[j]

        # Find which cell this point is in
        cell = searchsortedlast(xc, xj)

        # Determine local Riemann problem (between cell-1 and cell)
        if cell == 0
            # Left of domain: just take left state
            ρ_sol[j], u_sol[j], p_sol[j] = ρ[1], u[1], p[1]
            continue
        elseif cell == N
            # Right of domain
            ρ_sol[j], u_sol[j], p_sol[j] = ρ[N], u[N], p[N]
            continue
        end

        # Interface position
        xif = 0.5 * (xc[cell] + xc[cell+1])

        # Left and right states
        ρL, uL, pL = ρ[cell], u[cell], p[cell]
        ρR, uR, pR = ρ[cell+1], u[cell+1], p[cell+1]

        # Solve Riemann problem for this interface
        pstar, ustar = solve_star_region(ρL, uL, pL, ρR, uR, pR, γ)

        ξ = (xj - xif) / t
        ρj, uj, pj = sample_exact(ξ, ρL, uL, pL, ρR, uR, pR, pstar, ustar, γ)

        ρ_sol[j] = ρj
        u_sol[j] = uj
        p_sol[j] = pj
    end

    return ρ_sol, u_sol, p_sol
end