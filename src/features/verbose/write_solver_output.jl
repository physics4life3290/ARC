




function write_solver_output(log_file, W, U, F)
    write(log_file, "The output here is:

    Primitive Variables:
        Density Centers: $(W.density_centers)
        Velocity Centers: $(W.velocity_centers)
        Pressure Centers: $(W.pressure_centers)
        Internal Energy Centers: $(W.internal_energy_centers)

    Conserved Variables:
        Density Faces: $(U.density_centers)
        Momentum Faces: $(U.momentum_centers)
        Total Energy Faces: $(U.total_energy_centers)
    
    Flux Variables:
        Density Flux: $(F.density_flux)
        Momentum Flux: $(F.momentum_flux)
        Total Energy Flux: $(F.total_energy_flux)
    ")
    close(log_file)
end