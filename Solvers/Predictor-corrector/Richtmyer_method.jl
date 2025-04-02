


#   ∂U/∂t + ∂F/∂x = S 




function richtmyer_step(var_1, var_2, flux_var, grid, dx, dt, source_var)
    
    iters = length(var_1)
    new_var = zeros(iters)
    var_interp(var, x) = cubic_spline_interp(grid, var, x)
    
    if source_var == "none"
        for i in 2:iters-1
            r_interface = 0.5 * (var_1[i+1] + var_1[i]) - (dt/(4 * dx)) * (flux_var[i+1] - flux_var[i])
            l_interface = 0.5 * (var_1[i] + var_1[i-1]) - (dt/(4 * dx)) * (flux_var[i] - flux_var[i-1])

            r_var = var_interp(var_2, grid[i] + dx / 2)
            l_var = var_interp(var_2, grid[i] - dx / 2)
            new_var[i] = var_1[i] - (dt/dx) * (r_interface * r_var - l_interface * l_var)
        end
        return new_var
    elseif source_var != "none"
        for i in 2:iters-1
            r_interface = 0.5 * (var_1[i+1] + var_1[i]) - (dt/(4 * dx)) * (flux_var[i+1] - flux_var[i]) - (dt/(4 * dx)) * (source_var[i+1] - source_var[i])
            l_interface = 0.5 * (var_1[i] + var_1[i-1]) - (dt/(4 * dx)) * (flux_var[i] - flux_var[i-1]) - (dt/(4 * dx)) * (source_var[i] - source_var[i-1])
            
            r_var = var_interp(var_2, grid[i] + dx / 2)
            l_var = var_interp(var_2, grid[i] - dx / 2)
            new_var[i] = var_1[i] - (dt/dx) * (r_interface * r_var - l_interface * l_var)
        end
        return new_var
    end

end