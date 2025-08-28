





using Printf, Plots

# ---------- Sedov blast (spherical) ----------
γ   = 1.4
xL, xR = -1.0, 1.0
nx  = 800
x0  = 0.0
x   = range(xL, xR; length=nx)

# Ambient (pre-shock) state
ρ0, u0, p0 = 1.0, 0.0, 0.0          # cold medium

# Choose a target shock radius at t = tmax to fit in box
tmin = 1e-5
tmax = 1.0
Rtarget = 0.9 * (xR - xL)/2         # ~90% to boundary at tmax

# Sedov scaling (spherical): Rs(t) = C * t^(2/5)
# We don't actually need E or the exact similarity constant; we set C so Rs(tmax)=Rtarget
Csedov = Rtarget / (tmax^(2/5))

# Helper: shock radius and speed
Rs(t) = Csedov * t^(2/5)
Ss(t) = (2/5) * Csedov * t^(-3/5)   # d(Rs)/dt

# Strong-shock (Rankine–Hugoniot) jumps for γ-gas, upstream (ρ0,u0,p0) and shock speed S
# Immediate post-shock state (self-similar interior is more complex; we use front values here)
ρ2_of(S) = ((γ+1)/(γ-1)) * ρ0
p2_of(S) = 2ρ0 * S^2 / (γ+1)
u2_of(S) = 2S / (γ+1)               # gas velocity just behind the front, lab frame

# ---------- Animation ----------
nframes = 120
anim = @animate for k in 0:nframes
    t = tmin + (tmax - tmin) * (k / nframes)
    R = Rs(t)
    S = Ss(t)

    ρ = fill(ρ0, nx)
    u = fill(u0, nx)
    p = fill(p0, nx)

    # simple radial 1-D slice through a spherical blast, r = |x - x0|
    @inbounds for i in eachindex(x)
        r = abs(x[i] - x0)
        if r ≤ R
            # fill shocked region with immediate post-shock state
            ρ[i] = ρ2_of(S)
            u[i] = sign(x[i]-x0) * u2_of(S)  # outward on both sides along the line
            p[i] = p2_of(S)
        end
    end

    # plot three panels stacked (ρ, u, p)
    plot(size=(900,900), legend=false,
         left_margin=10Plots.mm, bottom_margin=10Plots.mm)

    plot!(x, ρ, xlabel="x", ylabel="ρ",
          title=@sprintf("Sedov blast (γ=%.1f), t=%.4f, R_s=%.3f", γ, t, R))
    vline!([-R + x0, R + x0], linestyle=:dash, alpha=0.6)

    plot!(x, u, xlabel="x", ylabel="u")
    vline!([-R + x0, R + x0], linestyle=:dash, alpha=0.6)

    plot!(x, p, xlabel="x", ylabel="p")
    vline!([-R + x0, R + x0], linestyle=:dash, alpha=0.6)
end

gif(anim, "sedov_strongshock.gif"; fps=30, show_msg=true)
println("Wrote sedov_strongshock.gif")
