






function lax_friedrichs_step_gpu(U::CuArray{T,2}, Δt::T, Δx::T, γ::T) where {T<:AbstractFloat}
    N = size(U, 2)
    Unew = similar(U)

    function lax_kernel!(U, Unew, Δt, Δx, γ, N)
        i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if 2 <= i <= N-1
            @inbounds begin
                # Manually extract components (no slicing!)
                ρL   = U[1, i-1]
                ρuL  = U[2, i-1]
                EL   = U[3, i-1]
                uL   = ρuL / ρL
                PL   = (γ - one(T)) * (EL - 0.5 * ρL * uL^2)
                FL1  = ρuL
                FL2  = ρuL * uL + PL
                FL3  = uL * (EL + PL)

                ρR   = U[1, i+1]
                ρuR  = U[2, i+1]
                ER   = U[3, i+1]
                uR   = ρuR / ρR
                PR   = (γ - one(T)) * (ER - 0.5 * ρR * uR^2)
                FR1  = ρuR
                FR2  = ρuR * uR + PR
                FR3  = uR * (ER + PR)

                Unew[1, i] = 0.5 * (ρR + ρL) - (Δt / (2Δx)) * (FR1 - FL1)
                Unew[2, i] = 0.5 * (ρuR + ρuL) - (Δt / (2Δx)) * (FR2 - FL2)
                Unew[3, i] = 0.5 * (ER + EL) - (Δt / (2Δx)) * (FR3 - FL3)
            end
        end
        return
    end

    threads = 256
    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks lax_kernel!(U, Unew, Δt, Δx, γ, N)

    # Apply boundary conditions on CPU
    Unew[:, 1] .= Unew[:, 2]
    Unew[:, N] .= Unew[:, N-1]

    return Unew
end