




function get_average(cell_interfaces, cell_centers, interface_variable, variable, interp_func)

    average_iters = length(variable)
    average = zeros(average_iters)

    for i in 1:average_iters
        integration_points = [cell_interfaces[i], cell_centers[i], cell_interfaces[i+1]]
        integrand_points = [interface_variable[i], variable[i], interface_variable[i+1]]
        integral = integrate(integration_points, integrand_points; interp_function=interp_func, weights=:simpsons, rule_type=:lownoise, stencil_size=6)
        average[i] = (1/(cell_interfaces[i+1] - cell_interfaces[i])) * integral
    end

    return average
end