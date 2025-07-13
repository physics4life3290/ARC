



struct PrimitiveVariables
    density::Vector{Float64}
    velocity::Vector{Float64}
    pressure::Vector{Float64}
    int_energy::Vector{Float64}
end

struct ConservativeVars
    density::Vector{Float64}
    momentum::Vector{Float64}
    total_energy::Vector{Float64}
end

function Construct1DConservatives(W::PrimitiveVariables, γ::Float64)
    U = ConservativeVars(
        W.density,
        W.density .* W.velocity,
        W.density .* (W.int_energy .+ W.pressure ./ (W.density .* (γ - 1)) .+ W.velocity.^2 ./ 2)
    )
    return U
end

struct FluxVars
    dens_flux::Vector{Float64}
    mome_flux::Vector{Float64}
    tot_ener_flux::Vector{Float64}
end