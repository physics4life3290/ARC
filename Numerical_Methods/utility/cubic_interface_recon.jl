




function cubic_reconstruction(cell_centers, cell_interfaces, avg_var)

    reconstruction_iters = length(cell_interfaces)
    recon_l_interface_var = zeros(reconstruction_iters)
    recon_r_interface_var = zeros(reconstruction_iters)


    for i in 3:reconstruction_iters-2
       
        l_center_points = [cell_centers[i-2], cell_centers[i-1], cell_centers[i]]
        l_center_values = [avg_var[i-2], avg_var[i-1], avg_var[i]]
        recon_l_interface_var[i] = cubic_spline_interp(l_center_points, l_center_values, cell_interfaces[i])

        r_center_points = [cell_centers[i-1], cell_centers[i], cell_centers[i+1]]
        r_center_values = [avg_var[i-1], avg_var[i], avg_var[i+1]]
        recon_r_interface_var[i] = cubic_spline_interp(r_center_points, r_center_values, cell_interfaces[i]) 
    end
    return recon_l_interface_var, recon_r_interface_var
end