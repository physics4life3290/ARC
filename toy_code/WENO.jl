using Plots

# --- Constants ---
const γ = 1.4
const nx = 100
const x = range(-1.0, 1.0, length=nx)
const dx = step(x)
const t_final = 5.0
const ε = 1e-36
const cfl = 0.3

# --- EOS ---
get_p(ρ, m, E) = (γ - 1.0) * (E - 0.5 * m^2 / ρ)

function get_flux(ρ, m, E)
    p = get_p(ρ, m, E)
    return [m, m^2/ρ + p, (E + p) * m / ρ]
end

# --- Characteristic Decomposition ---
function eigen_decomp(ρ, m, E)
    u = m / ρ
    p = get_p(ρ, m, E)
    a = sqrt(γ * p / ρ)
    H = (E + p) / ρ

    b1 = (γ - 1.0) / (2.0 * a^2)
    b2 = u / (2.0 * a)

    R = [1.0       1.0          1.0;
         u - a     u            u + a;
         H - u*a   0.5*u^2      H + u*a]

    L = [b1*0.5*u^2 + b2    -b1*u - 1/(2a)    b1;
         1.0 - b1*u^2        b1*u             -b1;
         b1*0.5*u^2 - b2    -b1*u + 1/(2a)    b1]

    return R, L, a
end

# --- WENO5 Left-Biased (i+1/2 from left) ---
@inline function weno5_left(v1, v2, v3, v4, v5)
    f1 =  (1/3)*v1 - (7/6)*v2 + (11/6)*v3
    f2 = -(1/6)*v2 + (5/6)*v3 +  (1/3)*v4
    f3 =  (1/3)*v3 + (5/6)*v4 -  (1/6)*v5

    β1 = (13/12)*(v1 - 2v2 + v3)^2 + (1/4)*(v1 - 4v2 + 3v3)^2
    β2 = (13/12)*(v2 - 2v3 + v4)^2 + (1/4)*(v2 - v4)^2
    β3 = (13/12)*(v3 - 2v4 + v5)^2 + (1/4)*(3v3 - 4v4 + v5)^2

    α1, α2, α3 = 0.1/(β1+ε)^2, 0.6/(β2+ε)^2, 0.3/(β3+ε)^2
    return (α1*f1 + α2*f2 + α3*f3) / (α1 + α2 + α3)
end

# --- WENO5 Right-Biased (i+1/2 from right) ---
@inline function weno5_right(v1, v2, v3, v4, v5)
    f1 = -(1/6)*v1 + (5/6)*v2 +  (1/3)*v3
    f2 =  (1/3)*v2 + (5/6)*v3 -  (1/6)*v4
    f3 = (11/6)*v3 - (7/6)*v4 +  (1/3)*v5

    β1 = (13/12)*(v1 - 2v2 + v3)^2 + (1/4)*(v1 - 4v2 + 3v3)^2
    β2 = (13/12)*(v2 - 2v3 + v4)^2 + (1/4)*(v2 - v4)^2
    β3 = (13/12)*(v3 - 2v4 + v5)^2 + (1/4)*(3v3 - 4v4 + v5)^2

    α1, α2, α3 = 0.1/(β1+ε)^2, 0.6/(β2+ε)^2, 0.3/(β3+ε)^2
    return (α1*f1 + α2*f2 + α3*f3) / (α1 + α2 + α3)
end

# --- Reflecting Ghost Cells ---

function apply_reflecting_ghosts(U)
    ρ, m, E = U[1], U[2], U[3]

    ρ_e = [ρ[5]; ρ[4]; ρ[3]; ρ[2]; ρ; ρ[end-1]; ρ[end-2]; ρ[end-3]; ρ[end-4]]
    m_e = [-m[5]; -m[4]; -m[3]; -m[2]; m; -m[end-1]; -m[end-2]; -m[end-3]; -m[end-4]]
    E_e = [E[5]; E[4]; E[3]; E[2]; E; E[end-1]; E[end-2]; E[end-3]; E[end-4]]

    return [ρ_e, m_e, E_e]
end

# --- RHS with Characteristic WENO5 ---
function compute_rhs(U)
    Ue = apply_reflecting_ghosts(U)
    ρe, me, Ee = Ue[1], Ue[2], Ue[3]

    # ghost offset: interior index i maps to extended index i+4
    ng = 4

    fluxes = [zeros(nx + 1) for _ in 1:3]

    @inbounds for i in 1:nx+1
        # Extended indices for left (L) and right (R) of interface i-1/2 in interior indexing
        # Interface between interior cells i-1 and i → extended index (i-1)+ng and i+ng
        ie = i + ng - 1  # left cell extended index (cell i-1 in 1-based interior)

        ρL, ρR = ρe[ie],   ρe[ie+1]
        mL, mR = me[ie],   me[ie+1]
        EL, ER = Ee[ie],   Ee[ie+1]

        # Roe-averaged state
        sqρL  = sqrt(ρL);  sqρR = sqrt(ρR)
        denom = sqρL + sqρR
        uL, uR = mL/ρL, mR/ρR
        HL = (EL + get_p(ρL, mL, EL)) / ρL
        HR = (ER + get_p(ρR, mR, ER)) / ρR

        u_roe = (sqρL*uL + sqρR*uR) / denom
        H_roe = (sqρL*HL + sqρR*HR) / denom
        a_roe = sqrt(abs((γ - 1.0) * (H_roe - 0.5*u_roe^2)))
        ρ_roe = sqρL * sqρR
        m_roe = ρ_roe * u_roe

        # FIX 1: Correct E_roe from Roe-averaged enthalpy H = (E+p)/ρ = γ/(γ-1)*p/ρ + u²/2
        # Inverting: p_roe/ρ_roe = (H_roe - 0.5*u_roe²) * (γ-1)/γ
        # E_roe = p_roe/(γ-1) + 0.5*ρ_roe*u_roe²
        p_roe = ρ_roe * (H_roe - 0.5*u_roe^2) * (γ - 1.0) / γ
        E_roe = p_roe / (γ - 1.0) + 0.5*m_roe^2 / ρ_roe

        # FIX 1 (cont): eigen_decomp now uses consistent Roe-averaged state, not EL
        R, L, _ = eigen_decomp(ρ_roe, m_roe, E_roe)

        # FIX 3: Compute λ_max over the full 6-point stencil, not just the interface pair.
        # This ensures the Lax-Friedrichs splitting is conservative for all stencil points.
        λ_max = 0.0
        for s in 1:6
            idx = ie - 3 + s  # stencil points: ie-2, ie-1, ie, ie+1, ie+2, ie+3
            ρs = ρe[idx]; ms = me[idx]; Es = Ee[idx]
            us = ms / ρs
            as = sqrt(γ * get_p(ρs, ms, Es) / ρs)
            λ_max = max(λ_max, abs(us) + as)
        end

        # Build characteristic flux splitting across 6-point stencil
        fp_char = [zeros(6) for _ in 1:3]
        fm_char = [zeros(6) for _ in 1:3]

        for s in 1:6
            idx = ie - 3 + s
            Ui  = [ρe[idx], me[idx], Ee[idx]]
            Fi  = get_flux(ρe[idx], me[idx], Ee[idx])

            fp_phys = 0.5 .* (Fi .+ λ_max .* Ui)
            fm_phys = 0.5 .* (Fi .- λ_max .* Ui)

            fp_c = L * fp_phys
            fm_c = L * fm_phys

            for k in 1:3
                fp_char[k][s] = fp_c[k]
                fm_char[k][s] = fm_c[k]
            end
        end

        # Reconstruct and project back
        flux_char = zeros(3)
        for k in 1:3
            fp_rec = weno5_left( fp_char[k][1], fp_char[k][2], fp_char[k][3],
                                  fp_char[k][4], fp_char[k][5])
            fm_rec = weno5_right(fm_char[k][2], fm_char[k][3], fm_char[k][4],
                                  fm_char[k][5], fm_char[k][6])
            flux_char[k] = fp_rec + fm_rec
        end

        flux_cons = R * flux_char

        for k in 1:3
            fluxes[k][i] = flux_cons[k]
        end
    end

    rhs = [zeros(nx) for _ in 1:3]
    @inbounds for k in 1:3
        for i in 1:nx
            rhs[k][i] = -(fluxes[k][i+1] - fluxes[k][i]) / dx
        end
    end

    return rhs
end

# --- SSP-RK3 ---
function ssp_rk3_step(U, dt)
    R1 = compute_rhs(U)
    U1 = [U[k] .+ dt .* R1[k] for k in 1:3]

    R2 = compute_rhs(U1)
    U2 = [0.75.*U[k] .+ 0.25.*U1[k] .+ 0.25.*dt.*R2[k] for k in 1:3]

    R3 = compute_rhs(U2)
    return [(1/3).*U[k] .+ (2/3).*U2[k] .+ (2/3).*dt.*R3[k] for k in 1:3]
end

# --- Positivity Check ---
function check_physical(U)
    any(U[1] .<= 0.0) && error("Negative density detected!")
    p = get_p.(U[1], U[2], U[3])
    any(p .<= 0.0) && error("Negative pressure detected!")
end

# --- Initial Conditions ---
function sod_ic()
    ρ = [xi < 0.0 ? 1.0  : 0.125 for xi in x]
    p = [xi < 0.0 ? 1.0  : 0.1   for xi in x]
    m = zeros(nx)
    E = p ./ (γ - 1.0) .+ 0.5 .* m.^2 ./ ρ
    return [ρ, m, E]
end

# --- Main Loop ---
function run_simulation(U; plot_every=0.05)
    current_t  = 0.0
    next_frame = 0.0
    frames     = []

    while current_t < t_final
        p  = get_p.(U[1], U[2], U[3])
        a  = sqrt.(γ .* p ./ U[1])
        u  = U[2] ./ U[1]
        dt = cfl * dx / maximum(abs.(u) .+ a)
        dt = min(dt, t_final - current_t)

        check_physical(U)

        if current_t >= next_frame
            push!(frames, (copy(U[1]), copy(p), copy(u), round(current_t, digits=3)))
            next_frame += plot_every
        end

        U = ssp_rk3_step(U, dt)
        current_t += dt
    end

    return U, frames
end

# --- Run ---
U0 = sod_ic()
U_final, frames = run_simulation(U0)

anim = @animate for (ρ_f, p_f, u_f, t) in frames
    plot(x, ρ_f, label="Density", lw=2,
         title="Characteristic WENO5 | Reflective Walls | t=$t",
         xlabel="x", ylims=(-0.5, 1.5))
    plot!(x, p_f, label="Pressure", lw=2)
    plot!(x, u_f, label="Velocity", lw=2)
end

mp4(anim, "weno5_characteristic.mp4", fps=20)