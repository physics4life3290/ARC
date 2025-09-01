




const MaybeVector{T} = Union{Nothing, Vector{T}}

struct PrimitiveVariables 

    density_centers::Union{Vector{Float64}, Float64, Nothing}
    velocity_centers::Union{Vector{Float64}, Float64, Nothing}
    pressure_centers::Union{Vector{Float64}, Float64, Nothing}
    internal_energy_centers::MaybeVector{Float64}
    density_faces::MaybeVector{Float64}
    velocity_faces::MaybeVector{Float64}
    pressure_faces::MaybeVector{Float64}
    internal_energy_faces::MaybeVector{Float64}
    
end

struct ConservativeVariables

    density_centers::Vector{Float64}
    momentum_centers::Vector{Float64}
    total_energy_centers::Vector{Float64}
    density_faces::MaybeVector{Float64}
    momentum_faces::MaybeVector{Float64}
    total_energy_faces::MaybeVector{Float64}

end

struct FluxVariables
    density_flux::Vector{Float64}
    momentum_flux::Vector{Float64}
    total_energy_flux::Vector{Float64}
end