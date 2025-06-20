using CUDA

function richtmyer_step_gpu!(U::CuDeviceMatrix, half_step_U::CuDeviceMatrix, Δt, Δx, γ)
    N = size(U, 2)
    threads = 256
    blocks = cld(N, threads)

    # Kernel 1: Compute half_step_U[:, i] for i=2:N-1
    function kernel_half_step!(U, half_step_U, Δt, Δx, N)
        i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if 2 ≤ i ≤ N-1
            u_im1 = @view U[:, i-1]
            u_ip1 = @view U[:, i+1]

            f_ip1 = flux_gpu(u_ip1, γ)
            f_im1 = flux_gpu(u_im1, γ)

            @inbounds for j in 1:size(U,1)
                half_step_U[j, i] = 0.5 * (u_ip1[j] + u_im1[j]) - (Δt / (2 * Δx)) * (f_ip1[j] - f_im1[j])
            end
        end
        return
    end

    # Kernel 2: Fix edges for half_step_U
    function kernel_fix_edges!(half_step_U, N)
        if threadIdx().x == 1 && blockIdx().x == 1
            @inbounds for j in 1:size(half_step_U, 1)
                half_step_U[j, 1] = half_step_U[j, 2]
                half_step_U[j, N] = half_step_U[j, N-1]
            end
        end
        return
    end

    # Kernel 3: Compute updated solution into half_step_U[:, i]
    function kernel_update_U!(U, half_step_U, Δt, Δx, N)
        i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if 2 ≤ i ≤ N-1
            u_i = @view U[:, i]
            hs_ip1 = @view half_step_U[:, i+1]
            hs_im1 = @view half_step_U[:, i-1]

            f_ip1 = flux_gpu(hs_ip1, γ)
            f_im1 = flux_gpu(hs_im1, γ)

            @inbounds for j in 1:size(U, 1)
                half_step_U[j, i] = u_i[j] - (Δt / Δx) * (f_ip1[j] - f_im1[j])
            end
        end
        return
    end

    # Kernel 4: Simple outflow BC for half_step_U
    function kernel_bc!(half_step_U, N)
        if threadIdx().x == 1 && blockIdx().x == 1
            @inbounds for j in 1:size(half_step_U, 1)
                half_step_U[j, 1] = half_step_U[j, 2]
                half_step_U[j, N] = half_step_U[j, N-1]
            end
        end
        return
    end

    # Launch kernels
    @cuda threads=threads blocks=blocks kernel_half_step!(U, half_step_U, Δt, Δx, N)
    @cuda threads=1 blocks=1 kernel_fix_edges!(half_step_U, N)
    @cuda threads=threads blocks=blocks kernel_update_U!(U, half_step_U, Δt, Δx, N)
    @cuda threads=1 blocks=1 kernel_bc!(half_step_U, N)

    # Copy the updated values back into U (overwrite in place)
    copyto!(U, half_step_U)

    return nothing
end
