




using Printf, Plots
include("../src/fluxes/riemann_solvers/ExactRiemannSolver.jl")
# ---------- Problem setup ----------
γ   = 1.4
xL, xR = -1.0, 1.0
nx  = 800
x0  = 0.0                 # initial discontinuity
x   = range(xL, xR; length=nx)

# Sod ICs
ρL, uL, pL = 1.0, 0.0, 1.0
ρR, uR, pR = 0.125, 0.0, 0.1

# Solve Riemann problem once (self-similar)
pstar, ustar = solve_star_pressure(ρL, uL, pL, ρR, uR, pR, γ)
ρLstar = star_density((ρL, uL, pL), pstar, γ)
ρRstar = star_density((ρR, uR, pR), pstar, γ)

println("Exact Riemann Solver:")
println("  p* = $(pstar)")
println("  u* = $(ustar)")
println("  ρ*_L = $(ρLstar)")
println("  ρ*_R = $(ρRstar)")

# ---------- Animation ----------
# We'll animate t in [tmin, tmax], sampling ξ = (x - x0)/t
tmin = 1e-6
tmax = 1.0
nframes = 120


anim = @animate for k in 0:nframes
    t = tmin + (tmax - tmin) * (k / nframes)
    ξ = [(xi - x0)/t for xi in x]

    ρ = similar(ξ)
    u = similar(ξ)
    p = similar(ξ)

    @inbounds for i in eachindex(ξ)
        ρi, ui, _pi = sample_exact(ξ[i], ρL, uL, pL, ρR, uR, pR, pstar, ustar, γ)
        ρ[i] = ρi
        u[i] = ui
        p[i] = _pi
    end
    plot(size=(900,900), xlabel="x", legend=false,
           left_margin=10Plots.mm, bottom_margin=10Plots.mm)
    plot!(x, ρ, xlabel="x", ylabel="ρ", title="Density ρ (t=$(round(t,digits=4)))")
    plot!(x, u, xlabel="x", ylabel="u", title="Velocity u")
    plot!(x, p, xlabel="x", ylabel="p", title="Pressure p")
    
end


gif(anim, "sod_exact.gif"; fps=30, show_msg=true)
println("Wrote sod_exact.gif")
