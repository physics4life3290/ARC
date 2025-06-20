using Plots
gr()

include("../src/build_conservative_vars.jl")
include("../src/fluxes.jl")
include("../src/time_step.jl")
include("../src/time_steppers/gpu/lax_friedrichs_gpu.jl")
include("../src/time_steppers/gpu/richtmyer_gpu.jl")
include("../src/riemann_solver/Riemann_Solver_GPU.jl")

function build_primitives(ρL, uL, PL, ρR, uR, PR, grid)
    W = zeros(Float64, 3, length(grid))
    for (j, x) in enumerate(grid)
        if x < 0
            W[:, j] .= (ρL, uL, PL)
        else
            W[:, j] .= (ρR, uR, PR)
        end
    end
    return W
end

function cons_to_prim(U::Matrix{Float64}, γ::Float64)
    W = zeros(size(U))  # Same shape: 3 x N
    for j in 1:size(U, 2)
        ρ = U[1, j]
        ρu = U[2, j]
        E = U[3, j]

        u = ρu / ρ
        p = (γ - 1.0) * (E - 0.5 * ρ * u^2)

        W[:, j] .= (ρ, u, p)
    end
    return W
end


function run_sod_shock_tube_gpu(ρL, uL, PL, ρR, uR, PR, grid, Δx, CFL, t_final)
    # === Initialization ===
    W  = build_primitives(ρL, uL, PL, ρR, uR, PR, grid)
    U = CuArray(Float32.(primitives_to_conservatives(W)))
    dt = compute_dt_gpu(U, Δx, CFL, γ)
    t  = 0.0f0
    step = 0

    # === Time Loop with Animation ===
    @gif for frame = 1:1000
        if t >= t_final
            break
        end

        #U = lax_friedrichs_step_gpu(U, dt, Δx, γ)
        #U = richtmyer_step_gpu!(U, dt, Δx, γ)
        #W,U = riemann_step!(W, grid, t, γ)
        t += dt
        dt = compute_dt_gpu(U, Δx, CFL, γ)
        step += 1

        U = Array(U)
        plt = plot(layout=(3,1), size=(800,900))  # Slightly larger size for safety

        plot!(plt[1], grid, U[1,:], label="", xlabel="x", ylabel="ρ", title="Density at t = $(round(t,digits=5)) s")
        plot!(plt[2], grid, U[2,:], label="", xlabel="x", ylabel="ρu", title="Momentum at t = $(round(t,digits=5)) s")
        plot!(plt[3], grid, U[3,:], label="", xlabel="x", ylabel="E", title="Energy at t = $(round(t,digits=5)) s")
        U = CuArray(U)  # Convert back to GPU array for next step
    end every 1
end




# === Parameters & Grid for cpu ===
const γ = Float32(4/3)
zones   = Int32(4096)
N       = zones
grid    = collect(range(-10.0f0, 10.0f0; length=N))  # Use Float32 for consistency
Δx      = (grid[end] - grid[1]) / N
CFL     = 0.3f0
t_final = 10.0f0

# === Initial Conditions ===
ρL, uL, PL = 1.0f0, 0.0f0, 1.0f0
ρR, uR, PR = 0.125f0, 0.0f0, 0.1f0



run_sod_shock_tube_gpu(ρL, uL, PL, ρR, uR, PR, grid, Δx, CFL, t_final)



#=
struct EvolutionScheme
    scheme::Symbol 
    gpu::Bool
    store_data::Symbol
end

struct sod_shock_tube{T}
    coord_sys::Union{Nothing, Symbol}
    zones::Int
    tube_length::T
    left_state::NTuple{3, T}
    right_state::NTuple{3, T}
    adiabatic_index::T
    solver::EvolutionScheme
    final_time::T
end

function SodShockTube(;
    coord_sys=nothing,
    zones::Int,
    tube_length,
    left_state,
    right_state,
    adiabatic_index,
    solver::EvolutionScheme,
    final_time
)
    T = typeof(tube_length)
    return sod_shock_tube{T}(
        coord_sys, zones, tube_length, left_state, right_state, adiabatic_index, solver, final_time
    )
end

solver = EvolutionScheme(
    :lax_friedrichs,
    false,
    :none
)

simulation = SodShockTube(
    coord_sys=nothing,
    zones=4096,
    tube_length=20.0,
    left_state=(1.0, 0.0, 1.0),
    right_state=(0.125, 0.0, 0.1),
    adiabatic_index=5/3,
    solver=solver,
    final_time=10.0
)
=#