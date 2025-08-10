function Construct1DBlastWavePrimitives(_grid, UserInput)

    dx = _grid.xcoord.spacing
    
    ρ = fill(ρ₀, nx)
    u = zeros(nx)
    p = fill(p₀, nx)

    # Center of the domain
    
    xc = parse(Float64, UserInput.secondary_input.grid_center)

    println(UserInput.secondary_input.blast_params)

    #=
    # Count cells inside the blast radius
    cells_in_blast = findall(abs.(x .- xc) .<= r₀)
    total_cells = length(cells_in_blast)

    # Energy per cell in the blast region
    ΔE = E₀ / total_cells

    # Convert to pressure using p = (γ - 1) * e, with e = energy density
    for i in cells_in_blast
        p[i] = (γ - 1.0) * (ΔE / dx)
    end

    W.density_centers .= ρ
    W.velocity_centers .= u
    W.pressure_centers .= p
    W.internal_energy_centers .= (p / ((UserInput.secondary_input.γ-1)*ρ))
    =#
end
