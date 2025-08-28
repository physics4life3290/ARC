



function apply_dirichlet_boundaries(ρ::AbstractVector, ρv::AbstractVector, E::AbstractVector, ng::Int, n::Int; ρf=nothing, ρvf=nothing, Ef=nothing)
    
    total = n + 2ng

    # Left boundary
    for i in 1:ng
        ρ[ng - i + 1]            = density_val
        ρv[ng - i + 1]           = momentum_val
        E[ng - i + 1]            = energy_val
    end

    # Right boundary
    for i in 1:ng
        ρ[total - ng + i]        = density_val
        ρv[total - ng + i]       = momentum_val
        E[total - ng + i]        = energy_val
    end

    if ρf !== nothing && ρvf !== nothing && Ef !== nothing
        total = length(ρf)

        for i in 1:ng
            ρf[ng - i + 1]       = density_val
            ρv[ng - i + 1]       = momentum_val
            Ef[ng - i + 1]       = energy_val
        end

        for i in 1:ng
            ρf[total - ng + i]   = density_val
            ρvf[total - ng + i]  = momentum_val
            Ef[total - ng + i]   = energy_val
        end
    elseif ρf === nothing || ρvf === nothing || Ef === nothing
        error("In 'src/grids/boundary_conditions/dirichlet.jl' all faces must be defined. Here is the current input-
        Density: $ρf
        Momentum: $ρvf
        Energy: $Ef")
    end

end
