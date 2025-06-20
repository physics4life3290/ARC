# exact_riemann.jl
# Standalone Julia script to solve and animate the exact Riemann problem for the 1D Euler equations.
# Requires Plots.jl for animation of density, velocity and pressure.

using Plots
gr()  # GR backend

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

# Solve star region
function solve_star_region(rhoL, uL, pL, rhoR, uR, pR, γ)
    p_old = guess_pressure(rhoL, uL, pL, rhoR, uR, pR, γ)
    for iter in 1:20
        fL, dFL = pressure_function(p_old, rhoL, pL, γ)
        fR, dFR = pressure_function(p_old, rhoR, pR, γ)
        f = fL + fR + (uR - uL)
        df = dFL + dFR
        p_new = p_old - f/df
        if abs(p_new - p_old) < 1e-6
            p_old = p_new
            break
        end
        p_old = max(1e-8, p_new)
    end
    p_star = p_old
    fL, _ = pressure_function(p_star, rhoL, pL, γ)
    fR, _ = pressure_function(p_star, rhoR, pR, γ)
    u_star = 0.5*(uL + uR + fR - fL)
    return p_star, u_star
end

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

# Exact Riemann solver
function exact_riemann(rhoL, uL, pL, rhoR, uR, pR; γ=1.4, x=0.0, t=1.0, nx=400)
    xi = LinRange(-10.0, 10.0, nx)
    sol = zeros(nx, 3)
    p_star, u_star = solve_star_region(rhoL, uL, pL, rhoR, uR, pR, γ)
    for i in 1:nx
        sol[i, :] .= sample_xi((xi[i] - x)/t, rhoL, uL, pL, rhoR, uR, pR, p_star, u_star, γ)
    end
    return xi, sol
end

function exact_riemann_profile(W::Matrix{Float64}, grid::Vector{Float64}; γ=1.4, t=1.0)
    nx = length(grid)
    sol = zeros(3, nx)  # variables × grid points, same shape as W

    # Extract left and right states from W (columns)
    ρL, uL, pL = W[:, 1]      # Leftmost state
    ρR, uR, pR = W[:, end]    # Rightmost state

    # Solve star region for single discontinuity at x=0
    p_star, u_star = solve_star_region(ρL, uL, pL, ρR, uR, pR, γ)

    # Loop over grid points
    for i in 1:nx
        ξ = grid[i] / t
        # sample_xi returns a 3-element vector [ρ, u, p] at ξ
        sol[:, i] .= sample_xi(ξ, ρL, uL, pL, ρR, uR, pR, p_star, u_star, γ)
    end

    return sol
end

function riemann_step!(W, grid, t)
    U = zeros(3, length(grid))
    W = exact_riemann_profile(W, Float64.(grid); t=t)
    U[1,:] = W[1,:]  # Update density
    U[2,:] = W[2,:] .* W[1,:]
    U[3,:] = W[3,:] ./ (γ - 1) .+ 0.5 .* W[1,:] .* W[2,:] .^ 2  # Update energy
    return W,U
end

#=
# Animation block: density, velocity, pressure
times = range(0.01, 10.0, length=1000)
anim = @animate for t in times
    xi, sol = exact_riemann(1.0, 0.0, 1.0, 0.125, 0.0, 0.1; γ=1.4, t=t, nx=400)
    plt1 = plot(xi, sol[:,1] ./ maximum(abs.(sol[:, 1])), ylim=(0,1.2), xlabel="x", ylabel="Density", title="t=$(round(t,digits=3))")
    plt2 = plot(xi, sol[:,2] ./ maximum(abs.(sol[:, 2])), ylim=(-1.0,1.0), xlabel="x", ylabel="Velocity")
    plt3 = plot(xi, sol[:,3] ./ maximum(abs.(sol[:, 3])), ylim=(0,1.5), xlabel="x", ylabel="Pressure")
    plot(plt1, plt2, plt3, layout = (3,1), link=:x)
end

# Save to GIF
gif(anim, "riemann.gif", fps=15)

# To run: julia exact_riemann.jl
=#