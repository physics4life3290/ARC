




@inline l2norm(v) = sqrt(sum(abs2, v))

function count_modified!(orig, newL, newR)
    c = 0
    @inbounds for i in eachindex(orig)
        if newL[i] != orig[i] || newR[i] != orig[i]
            c += 1
        end
    end
    return c
end

function GodunovStepDebug!(
    W::PrimitiveVariables, U::ConservativeVariables, F::FluxVariables,
    reconstruction::Symbol, limiter::Symbol,
    flattening::Bool, steepening::Bool,
    boundary_condition::Symbol, riemanntype::Symbol, 
    γ::Float64, spacing, dt::Float64,
    total_zones::Int32, zones::Int32, ghost_zones::Int32, grid_points::Vector{Float64},
    mode::Symbol,
    step::Int, t::Float64, CFL::Float64
)

    step_start = time()

    io = open("GodunovLog.txt", "a")

    # =========================
    # INITIAL NORMS
    # =========================
    ρ0 = l2norm(W.centers.density)
    mom0 = l2norm(U.centers.momentum[1])
    E0 = l2norm(U.centers.total_energy)

    recon_time = 0.0
    limit_time = 0.0
    riemann_time = 0.0
    update_time = 0.0

    flattened_cells = 0
    steepened_cells = 0

    max_dp = 0.0
    max_wave = 0.0
    shock_count = 0

    # =========================
    # RECONSTRUCTION
    # =========================
    recon_time = @elapsed begin

        if reconstruction != :Constant
            if reconstruction == :Linear
                slope_ρ = compute_slopes(W.centers.density, limiter)
                slope_u = compute_slopes(W.centers.velocity[1], limiter)
                slope_p = compute_slopes(W.centers.pressure, limiter)

                ρL, ρR = linear_reconstruct(W.centers.density, slope_ρ)
                uL, uR = linear_reconstruct(W.centers.velocity[1], slope_u)
                pL, pR = linear_reconstruct(W.centers.pressure, slope_p)

            elseif reconstruction == :Parabolic
                ρL, ρR = parabolic_reconstruct(W.centers.density)
                uL, uR = parabolic_reconstruct(W.centers.velocity[1])
                pL, pR = parabolic_reconstruct(W.centers.pressure)

            elseif reconstruction == :Cubic
                ρL, ρR = reconstruct_interfaces(W.centers.density, grid_points)
                uL, uR = reconstruct_interfaces(W.centers.velocity[1], grid_points)
                pL, pR = reconstruct_interfaces(W.centers.pressure, grid_points)
            end
        else
            ρL, ρR = W.centers.density, W.centers.density
            uL, uR = W.centers.velocity[1], W.centers.velocity[1]
            pL, pR = W.centers.pressure, W.centers.pressure
        end
    end

    # =========================
    # LIMIT / FLATTEN / STEEPEN
    # =========================
    limit_time = @elapsed begin

        if steepening
            steepened_cells = count_modified!(W.centers.density, ρL, ρR)
            ρL, ρR = steepening_contact!(
                W.centers.density, W.centers.pressure,
                grid_points, γ, ρL, ρR
            )
        end

        if flattening
            flat_coeff = compute_flattening_coeff!(
                zeros(length(W.centers.density)),
                W.centers.pressure, W.centers.velocity[1]
            )

            flattened_cells = count(x -> x < 0.999, flat_coeff)

            @inbounds begin
                ρL .= W.centers.density .+ flat_coeff .* (ρL .- W.centers.density)
                ρR .= W.centers.density .+ flat_coeff .* (ρR .- W.centers.density)
                uL .= W.centers.velocity[1] .+ flat_coeff .* (uL .- W.centers.velocity[1])
                uR .= W.centers.velocity[1] .+ flat_coeff .* (uR .- W.centers.velocity[1])
                pL .= W.centers.pressure .+ flat_coeff .* (pL .- W.centers.pressure)
                pR .= W.centers.pressure .+ flat_coeff .* (pR .- W.centers.pressure)
            end
        end
    end

    # =========================
    # RIEMANN + FLUXES
    # =========================
    riemann_time = @elapsed begin
        for i in 2:total_zones
            @inbounds begin
                f1 = f2 = f3 = 0.0

                dp = abs(pL[i-1] - pR[i])
                max_dp = max(max_dp, dp)
                if dp > 1e-3
                    shock_count += 1
                end

                if riemanntype == :HLL
                    f1, f2, f3 = Riemann_HLL(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)

                elseif riemanntype == :HLLC
                    f1, f2, f3 = Riemann_HLLC(ρL[i-1], uL[i-1], pL[i-1], ρR[i], uR[i], pR[i], γ)

                elseif riemanntype == :Exact
                    pstar, ustar = solve_star_pressure(
                        ρL[i-1], uL[i-1], pL[i-1],
                        ρR[i], uR[i], pR[i], γ
                    )
                    prims = sample_exact(0.0,
                        ρL[i-1], uL[i-1], pL[i-1],
                        ρR[i], uR[i], pR[i],
                        pstar, ustar, γ
                    )

                    f1 = prims[1] * prims[2]
                    f2 = prims[1] * prims[2]^2 + prims[3]
                    E = prims[3]/(γ-1) + 0.5*prims[1]*prims[2]^2
                    f3 = prims[2] * (E + prims[3])

                    max_wave = max(max_wave, abs(ustar))
                end

                F.density_flux[i] = f1
                F.momentum_flux[1][i] = f2
                F.total_energy_flux[i] = f3
            end
        end
    end

    # =========================
    # UPDATE STEP
    # =========================
    update_time = @elapsed begin
        for i in 2:total_zones-1
            @inbounds begin
                U.centers.density[i] -= dt/spacing * (F.density_flux[i+1] - F.density_flux[i])
                U.centers.momentum[1][i] -= dt/spacing * (F.momentum_flux[1][i+1] - F.momentum_flux[1][i])
                U.centers.total_energy[i] -= dt/spacing * (F.total_energy_flux[i+1] - F.total_energy_flux[i])
            end
        end
    end

    # =========================
    # FINAL NORMS
    # =========================
    ρf = l2norm(U.centers.density)
    momf = l2norm(U.centers.momentum[1])
    Ef = l2norm(U.centers.total_energy)

    step_time = time() - step_start

    # =========================
    # WRITE LOG
    # =========================
    println(io, "Step: $step | t = $t | Δt = $dt | CFL = $CFL | Step time = $(round(step_time, digits=4)) s\n")

    println(io, "Variable: Density (ρ)")
    println(io, "  Norms: initial=$ρ0 | post-step=$ρf")
    println(io, "  Limiter: $limiter | Flattened: $flattened_cells cells | Steepened: $steepened_cells cells\n")

    println(io, "Variable: Momentum (ρu)")
    println(io, "  Norms: initial=$mom0 | post-step=$momf\n")

    println(io, "Variable: Total Energy (E)")
    println(io, "  Norms: initial=$E0 | post-step=$Ef\n")

    println(io, "Timing (s): reconstruction=$recon_time | flatten/steepen=$limit_time")
    println(io, "Riemann diagnostics: max Δp = $max_dp | shock count = $shock_count | max wave speed = $max_wave | timing=$riemann_time s")
    println(io, "Step update timing = $update_time s\n")

    close(io)
end