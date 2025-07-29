






function ExactRiemannSolve(states, point, t; γ=1.4)
    ρL, uL, pL, ρR, uR, pR = states
    # Solve star region for single discontinuity at x=0
    p_star, u_star = solve_star_region(ρL, uL, pL, ρR, uR, pR, γ)
    ξ = point / t
    ρ, u, p = sample_xi(ξ, ρL, uL, pL, ρR, uR, pR, p_star, u_star, γ)

    return ρ, u, p
end