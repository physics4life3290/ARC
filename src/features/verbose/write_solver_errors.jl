



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