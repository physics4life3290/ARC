




function flux(Ui::AbstractVector)
    ρ, ρu, E = Ui
    u = ρu / ρ
    P = (γ - 1)*(E - 0.5*ρ*u^2)
    return [ ρu,
             ρu*u + P,
             u*(E + P) ]
end

# Flux kernel (same as before)
@inline function flux_gpu(U::CuArray{Float32, 2}, γ::Float32)
    N = size(U, 2)
    F = similar(U)
    threads = 256
    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks flux_kernel!(F, U, γ, N)
    return F
end

function flux_kernel!(
    F::CuDeviceMatrix{Float32},
    U::CuDeviceMatrix{Float32},
    γ::Float32,
    N::Int,
)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x + 1
    if i >= 1 && i <= N
        ρ = max(U[1, i], 1e-6f0)
        u = U[2, i] / ρ
        E = U[3, i]
        p = (γ - 1f0) * (E - 0.5f0 * ρ * u * u)

        F[1, i] = ρ * u
        F[2, i] = ρ * u^2 + p
        F[3, i] = u * (E + p)
    end
    return
end