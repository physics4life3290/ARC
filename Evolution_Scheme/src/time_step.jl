




function compute_dt(U, Δx, CFL)
    max_speed = 0.0
    for j in 1:size(U,2)
        ρ, ρu, E = U[:,j]
        u = ρu / ρ
        P = (γ - 1) * (E - 0.5 * ρ * u^2)
        c = sqrt(γ * abs(P / ρ))
        speed = abs(u) + c
        max_speed = max(max_speed, speed)
    end
    return CFL * Δx / max_speed
end

function compute_dt_gpu(U::CuArray{T,2}, Δx::T, CFL::T, γ::T) where {T<:AbstractFloat}
    N = size(U, 2)
    speeds = CUDA.zeros(T, N)

    # Kernel to compute wavespeed per zone
    function speed_kernel(U, speeds, γ, N)
        j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if j <= N
            @inbounds begin
                ρ  = U[1, j]
                ρu = U[2, j]
                E  = U[3, j]

                u = ρu / ρ
                u2 = u * u
                P = (γ - one(T)) * (E - 0.5 * ρ * u2)
                c = sqrt(γ * abs(P / ρ))
                speed = abs(u) + c
                speeds[j] = speed
            end
        end
        return
    end

    threads = 256
    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks speed_kernel(U, speeds, γ, N)

    # Reduce to get maximum speed
    max_speed = CUDA.maximum(speeds)

    return CFL * Δx / max_speed
end