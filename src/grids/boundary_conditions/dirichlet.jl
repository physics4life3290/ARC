




function apply_dirichlet_boundaries!(ρ::AbstractVector, ρv::AbstractVector, E::AbstractVector, ng::Int, n::Int;
                                    density_val::Float64=1.0,
                                    momentum_val::Float64=0.0,
                                    energy_val::Float64=1.0,
                                    ρf=nothing, ρvf=nothing, Ef=nothing)

    total = n + 2ng

    # Left boundary (centers)
    @inbounds @simd for i in 1:ng
        dest = ng - i + 1
        ρ[dest]  = density_val
        ρv[dest] = momentum_val
        E[dest]  = energy_val
    end

    # Right boundary (centers)
    @inbounds @simd for i in 1:ng
        dest = total - ng + i
        ρ[dest]  = density_val
        ρv[dest] = momentum_val
        E[dest]  = energy_val
    end

    # Faces (if provided)
    if ρf !== nothing && ρvf !== nothing && Ef !== nothing
        total_faces = length(ρf)

        # Left boundary (faces)
        @inbounds @simd for i in 1:ng
            dest = ng - i + 1
            ρf[dest]   = density_val
            ρvf[dest]  = momentum_val
            Ef[dest]   = energy_val
        end

        # Right boundary (faces)
        @inbounds @simd for i in 1:ng
            dest = total_faces - ng + i
            ρf[dest]   = density_val
            ρvf[dest]  = momentum_val
            Ef[dest]   = energy_val
        end
    elseif ρf === nothing || ρvf === nothing || Ef === nothing
        error("In 'apply_dirichlet_boundaries!' all faces must be defined. Got: ρf=$ρf, ρvf=$ρvf, Ef=$Ef")
    end

    return nothing
end

