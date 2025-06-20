




using CUDA
using Plots
gr()

# CPU functions (unchanged)
function guess_pressure(rhoL, uL, pL, rhoR, uR, pR, γ)
    cL = sqrt(γ * pL / rhoL)
    cR = sqrt(γ * pR / rhoR)
    p_pv = 0.5*(pL + pR) - 0.125*(uR - uL)*(rhoL + rhoR)*(cL + cR)
    return max(1e-6, p_pv)
end

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

# GPU kernel to sample all xi values
function sample_kernel!(sol, xi, rhoL, uL, pL, rhoR, uR, pR, p_star, u_star, γ)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    N = length(xi)
    if i > N
        return
    end

    x = xi[i]
    sρ, su, sp = 0.0, 0.0, 0.0
    xi_val = x

    if xi_val <= u_star
        if p_star > pL
            sL = uL - sqrt((γ+1)/(2γ)*p_star/pL + (γ-1)/(2γ))*sqrt(γ*pL/rhoL)
            if xi_val <= sL
                sρ, su, sp = rhoL, uL, pL
            else
                sρ = rhoL*((p_star/pL + (γ-1)/(γ+1))/((γ-1)/(γ+1)*p_star/pL + 1))
                su = u_star
                sp = p_star
            end
        else
            aL = sqrt(γ*pL/rhoL)
            shL = uL - aL
            stL = u_star - sqrt(γ*p_star/rhoL)
            if xi_val <= shL
                sρ, su, sp = rhoL, uL, pL
            elseif xi_val <= stL
                u = (2/(γ+1))*(aL + (γ-1)/2*uL + xi_val)
                a = (2/(γ+1))*(aL + (γ-1)/2*(uL - xi_val))
                sρ = rhoL*(a/aL)^(2/(γ-1))
                su = u
                sp = pL*(a/aL)^(2*γ/(γ-1))
            else
                sρ = rhoL*(p_star/pL)^(1/γ)
                su = u_star
                sp = p_star
            end
        end
    else
        if p_star > pR
            sR = uR + sqrt((γ+1)/(2γ)*p_star/pR + (γ-1)/(2γ))*sqrt(γ*pR/rhoR)
            if xi_val >= sR
                sρ, su, sp = rhoR, uR, pR
            else
                sρ = rhoR*((p_star/pR + (γ-1)/(γ+1))/((γ-1)/(γ+1)*p_star/pR + 1))
                su = u_star
                sp = p_star
            end
        else
            aR = sqrt(γ*pR/rhoR)
            shR = uR + aR
            stR = u_star + sqrt(γ*p_star/rhoR)
            if xi_val >= shR
                sρ, su, sp = rhoR, uR, pR
            elseif xi_val >= stR
                u = (2/(γ+1))*(-aR + (γ-1)/2*uR + xi_val)
                a = (2/(γ+1))*(aR - (γ-1)/2*(uR - xi_val))
                sρ = rhoR*(a/aR)^(2/(γ-1))
                su = u
                sp = pR*(a/aR)^(2*γ/(γ-1))
            else
                sρ = rhoR*(p_star/pR)^(1/γ)
                su = u_star
                sp = p_star
            end
        end
    end

    sol[i, 1] = sρ
    sol[i, 2] = su
    sol[i, 3] = sp
    return
end

# GPU-enabled Riemann solver
function exact_riemann_gpu(rhoL, uL, pL, rhoR, uR, pR; γ=1.4, x=0.0, t=1.0, nx=400)
    xi_cpu = LinRange(-10.0, 10.0, nx)
    xi = CUDA.fill(0.0f64, nx)
    @. xi = (xi_cpu - x) / t

    sol = CUDA.zeros(Float64, nx, 3)
    p_star, u_star = solve_star_region(rhoL, uL, pL, rhoR, uR, pR, γ)

    threads = 256
    blocks = cld(nx, threads)
    @cuda threads=threads blocks=blocks sample_kernel!(sol, xi, rhoL, uL, pL, rhoR, uR, pR, p_star, u_star, γ)

    return xi_cpu, Array(sol)
end

using CUDA

function exact_riemann_profile_GPU(W::CuArray{Float32,2}, grid::Vector{Float32}; γ::Float32=1.4f0, t::Float32=1.0f0)

    nx = length(grid)

    # Extract left and right states from W (columns)
    ρL, uL, pL = W[:, 1]
    ρR, uR, pR = W[:, end]

    # Solve star region (on CPU)
    p_star, u_star = solve_star_region(ρL, uL, pL, ρR, uR, pR, γ)

    # Prepare device arrays
    ξ_host = grid ./ t
    ξ_dev = CuArray(ξ_host)
    sol_dev = CUDA.zeros(Float64, nx, 3)

    # Launch GPU kernel
    threads = 256
    blocks = cld(nx, threads)
    @cuda threads=threads blocks=blocks sample_kernel!(
        sol_dev, ξ_dev, ρL, uL, pL, ρR, uR, pR, p_star, u_star, γ
    )

    # Transfer solution back to CPU
    sol = Array(sol_dev)'
    return sol  # shape (3, nx)
end

function riemann_step_gpu(W, grid, t; γ=1.4)
    W = exact_riemann_profile_GPU(W, Float64.(grid); t=t, γ=γ)

    U = similar(W)  # or define externally if needed
    U[1, :] = W[1, :]
    U[2, :] = W[2, :] .* W[1, :]
    U[3, :] = W[3, :] ./ (γ - 1) .+ 0.5 .* W[1, :] .* W[2, :].^2
    return W, U
end


#=
# Animation block
times = range(0.01, 10.0, length=300)
anim = @animate for t in times
    xi, sol = exact_riemann_gpu(1.0, 0.0, 1.0, 0.125, 0.0, 0.1; γ=1.4, t=t, nx=800)
    plt1 = plot(xi, sol[:,1] ./ maximum(abs.(sol[:, 1])), ylim=(0,1.2), xlabel="x", ylabel="Density", title="t=$(round(t,digits=3))")
    plt2 = plot(xi, sol[:,2] ./ maximum(abs.(sol[:, 2])), ylim=(-1.0,1.0), xlabel="x", ylabel="Velocity")
    plt3 = plot(xi, sol[:,3] ./ maximum(abs.(sol[:, 3])), ylim=(0,1.5), xlabel="x", ylabel="Pressure")
    plot(plt1, plt2, plt3, layout = (3,1), link=:x)
end

gif(anim, "riemann_gpu.gif", fps=15)
=#