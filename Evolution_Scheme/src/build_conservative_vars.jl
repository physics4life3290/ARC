


using CUDA

function primitives_to_conservatives(W)
    U = zeros(size(W))
    for j in 1:size(W,2)
        ρ, u, P = W[:,j]
        U[1,j] = ρ
        U[2,j] = ρ .* u
        U[3,j] = P ./ (γ - 1) .+ 0.5 .* ρ .* u .^2
    end
    return U
end

function primitives_to_conservatives_II(W)
    if ndims(W) == 1 || size(W,2) == 1
        ρ, u, P = W[:,1]
        return [
            ρ;
            ρ * u;
            P / (γ - 1) + 0.5 * ρ * u^2
        ]
    end

    U = zeros(size(W))
    for j in 1:size(W, 2)
        ρ, u, P = W[:, j]
        U[1, j] = ρ
        U[2, j] = ρ * u
        U[3, j] = P / (γ - 1) + 0.5 * ρ * u^2
    end
    return U
end


# GPU version
function primitives_to_conservatives_gpu(W::CuArray{T,2}, γ::T) where {T}
    U = similar(W)
    N = size(W, 2)

    function kernel(W, U, γ, N)
        j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if j <= N
            ρ = W[1, j]
            u = W[2, j]
            P = W[3, j]
            U[1, j] = ρ
            U[2, j] = ρ * u
            U[3, j] = P / (γ - 1) + 0.5f0 * ρ * u * u
        end
        return
    end

    threads = 256
    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks kernel(W, U, γ, N)
    return U
end