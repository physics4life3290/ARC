



# --- Total Variation (TV) ---
function total_variation(u)
    tv = 0.0
    for i in 1:(length(u)-1)
        tv += abs(u[i+1] - u[i])
    end
    return tv
end

function solver_error_metrics(var_old, var_new, var_flux, CFL, Δx, Δt, logfile)
    
    N = length(var_old)

    # --- Conservation check ---
    Q_old = sum(var_old)
    Q_new = sum(var_new)
    conservation_error = Q_new - Q_old


    TV_old = total_variation(var_old)
    TV_new = total_variation(var_new)
    TV_change = TV_new - TV_old

    # --- Residual (update norm) ---
    diff = var_new .- var_old
    residual_norm = sqrt(sum(diff.^2) / N)

    # --- Flux imbalance ---
    flux_in = var_flux[1]
    flux_out = var_flux[end]
    flux_imbalance = flux_in - flux_out

    # --- Logging ---
    
    ts = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    println(logfile, "\n========== Solver Error Diagnostics ==========")
    println(logfile, "Timestamp: $ts")
    println(logfile, @sprintf("Grid size: %d", N))
    println(logfile, @sprintf("Δx: %.6e, Δt: %.6e", Δx, Δt))
    println(logfile, @sprintf("CFL number: %.6e", CFL))
    println(logfile, @sprintf("Conservation error: %.6e", conservation_error))
    println(logfile, @sprintf("Total variation change: %.6e", TV_change))
    println(logfile, @sprintf("Residual norm: %.6e", residual_norm))
    println(logfile, @sprintf("Flux imbalance: %.6e", flux_imbalance))
    println(logfile, "=============================================\n")
    
    
end

function solver_diagnostics(U, U_old, F, CFL, Δx, Δt; logfile=nothing)
    N = length(U.density_centers)
    ts = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    # --- Conservation checks (mass example; can add momentum/energy if needed)
    ρ_old_sum = sum(U_old.density_centers)
    ρ_new_sum = sum(U.density_centers)
    conservation_error = ρ_new_sum - ρ_old_sum

    # --- Flux imbalance (per variable) ---
    flux_imbalance_ρ = F.density_flux[1]      - F.density_flux[end]
    flux_imbalance_m  = F.momentum_flux[1]    - F.momentum_flux[end]
    flux_imbalance_E  = F.total_energy_flux[1] - F.total_energy_flux[end]

    # --- Step increments ---
    Δρ = U.density_centers       .- U_old.density_centers
    Δm = U.momentum_centers      .- U_old.momentum_centers
    ΔE = U.total_energy_centers  .- U_old.total_energy_centers

    # --- Norms ---
    L2_ρ = sqrt(sum(Δρ .^ 2))
    L2_m = sqrt(sum(Δm .^ 2))
    L2_E = sqrt(sum(ΔE .^ 2))

    Linf_ρ, Linf_m, Linf_E = maximum(abs.(Δρ)), maximum(abs.(Δm)), maximum(abs.(ΔE))

    # --- Total variation ---
    TV_ρ_new = sum(abs.(diff(U.density_centers)))
    TV_m_new = sum(abs.(diff(U.momentum_centers)))
    TV_E_new = sum(abs.(diff(U.total_energy_centers)))

    TV_ρ_old = sum(abs.(diff(U_old.density_centers)))
    TV_m_old = sum(abs.(diff(U_old.momentum_centers)))
    TV_E_old = sum(abs.(diff(U_old.total_energy_centers)))

    TV_change = (TV_ρ_new + TV_m_new + TV_E_new) - (TV_ρ_old + TV_m_old + TV_E_old)

    # --- Max/Min values ---
    max_ρ, min_ρ = maximum(U.density_centers), minimum(U.density_centers)
    max_m, min_m = maximum(U.momentum_centers), minimum(U.momentum_centers)
    max_E, min_E = maximum(U.total_energy_centers), minimum(U.total_energy_centers)

    # --- Residual norm (global update) ---
    diff_all = vcat(Δρ, Δm, ΔE)
    residual_norm = sqrt(sum(diff_all.^2) / length(diff_all))

    # --- Print/log ---
    output = """
    ========== Solver Diagnostics ==========
    Timestamp: $ts
    Grid size: $N
    Δx: $(Δx), Δt: $(Δt)
    CFL number: $(CFL)

    -- Conservation --
    Mass conservation error: $conservation_error

    -- Flux Imbalance --
    Density flux imbalance: $flux_imbalance_ρ
    Momentum flux imbalance: $flux_imbalance_m
    Energy flux imbalance: $flux_imbalance_E

    -- L2 Norms (Δ) --
    Density: $L2_ρ
    Momentum: $L2_m
    Energy: $L2_E

    -- L∞ Norms (Δ) --
    Density: $Linf_ρ
    Momentum: $Linf_m
    Energy: $Linf_E

    -- Total Variation --
    TV (ρ,m,E): $TV_ρ_new, $TV_m_new, $TV_E_new
    Total variation change: $TV_change

    -- Max/Min Values --
    Density: ($max_ρ, $min_ρ)
    Momentum: ($max_m, $min_m)
    Energy: ($max_E, $min_E)

    -- Residual Norm --
    $residual_norm
    ========================================
    """

    if logfile !== nothing
        println(logfile, output)
    end
    
end

