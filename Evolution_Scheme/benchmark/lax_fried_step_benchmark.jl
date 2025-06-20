




println("Benchmarking Lax-Friedrichs method... ")
# Problem setup
ρL, uL, PL = 1.0f0, 0.0f0, 1.0f0
ρR, uR, PR = 0.125f0, 0.0f0, 0.1f0

W_cpu = build_primitives(ρL, uL, PL, ρR, uR, PR, grid)
W_gpu = CuArray(W_cpu)  # upload after building on CPU
println("Building Conservative Variables on the cpu...")
U_cpu = @btime primitives_to_conservatives(W_cpu)
println("Building Conservative Variables on the gpu...")
U_gpu = @btime primitives_to_conservatives_gpu(W_gpu, Float32(γ))

println("Starting Flux computation on CPU...")
F_cpu = zeros(size(U_cpu))
@btime for j in 1:size(U_cpu, 2)
    F_cpu[:, j] = flux(U_cpu[:, j])
end
println("Starting Flux computation on GPU...")
F_gpu = @btime flux_gpu(U_gpu, Float32(γ))
println("Flux computation completed.")

println("Computing time-step on CPU...")
dt_cpu = @btime compute_dt(U_cpu, Δx, CFL)
println("Time-step on CPU: $dt_cpu")
println("Computing time-step on GPU...")
dt_gpu = @btime compute_dt_gpu(U_gpu, Float32(Δx), Float32(CFL), Float32(γ))
println("Time-step on GPU: $dt_gpu")

println("Performing Lax-Friedrichs step on CPU...")
New_U_cpu = @btime lax_friedrichs_step!(U_cpu, dt_cpu, Δx)
println("Performing Lax-Friedrichs step on GPU...")
New_U_gpu = @btime lax_friedrichs_step_gpu(U_gpu, Float32(dt_gpu), Float32(Δx), Float32(γ))
println("Lax-Friedrichs step completed.")
using LinearAlgebra
norm_diff = norm(New_U_cpu .- Array(New_U_gpu))
println("Norm of difference between CPU and GPU results: $norm_diff")
