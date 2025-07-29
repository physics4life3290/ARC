




struct PrimitiveVariablesFVM
    density_centers::Vector{Float64}
    velocity_centers::Vector{Float64}
    pressure_centers::Vector{Float64}
    int_energy_centers::Vector{Float64}
    density_faces::Vector{Float64}
    velocity_faces::Vector{Float64}
    pressure_faces::Vector{Float64}
    int_energy_faces::Vector{Float64}
end

struct ConservativeVarsFVM
    density_centers::Vector{Float64}
    momentum_centers::Vector{Float64}
    total_energy_centers::Vector{Float64}
    density_faces::Vector{Float64}
    momentum_faces::Vector{Float64}
    total_energy_faces::Vector{Float64}
end

function Construct1DFVMConservatives(W::PrimitiveVariables, γ::Float64)
    U = ConservativeVarsFVM(
        W.density_centers,
        W.density_centers .* W.velocity_centers,
        W.density_centers .* (W.int_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (γ - 1)) .+ W.velocity_centers.^2 ./ 2),
        W.density_faces,
        W.density_faces .* W.velocity_faces,
        W.density_faces .* (W.int_energy_faces .+ W.pressure_faces ./ (W.density_faces .* (γ - 1)) .+ W.velocity_faces.^2 ./ 2)
    )
    return U
end

struct FluxVarsFVM
    dens_flux::Vector{Float64}
    mome_flux::Vector{Float64}
    tot_ener_flux::Vector{Float64}
end