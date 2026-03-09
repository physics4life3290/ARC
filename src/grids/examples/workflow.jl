struct PrimitiveFields{T,N}
    density::Vector{T}
    velocity::NTuple{N, Vector{T}}
    pressure::Vector{T}
end

struct PrimitiveVariables{T,N}
    centers::PrimitiveFields{T,N}
    faces::PrimitiveFields{T,N}
end

struct ConservativeFields{T,N}
    density::Vector{T}
    momentum::NTuple{N, Vector{T}}
    total_energy::Vector{T}
end

struct ConservativeVariables{T,N}
    centers::ConservativeFields{T,N}
    faces::ConservativeFields{T,N}
end

struct FluxVariables{T, N}
    density_flux::Vector{T}
    momentum_flux::NTuple{N, Vector{T}}
    total_energy_flux::Vector{T}
end

# ── EOS Interface ──────────────────────────────────────────────────────────────

abstract type AbstractEOS end

struct IdealGasEOS <: AbstractEOS
    γ::Float64
end

pressure(eos::IdealGasEOS, ρ, e)        = (eos.γ - 1) * ρ * e
sound_speed(eos::IdealGasEOS, ρ, e)     = sqrt(eos.γ * (eos.γ - 1) * e)
specific_energy(eos::IdealGasEOS, ρ, P) = P / ((eos.γ - 1) * ρ)
enthalpy(eos::IdealGasEOS, ρ, P, e)     = e + P / ρ

# Generic total_energy works for any EOS that implements specific_energy
function total_energy(eos::AbstractEOS, ρ, v, P)
    e = specific_energy(eos, ρ, P)
    return ρ * e + 0.5 * ρ * v^2
end

# Generic primitive recovery — works for any EOS that implements pressure
function prim_from_cons(eos::AbstractEOS, ρ, ρv, E)
    v = ρv / ρ
    e = (E - 0.5 * ρ * v^2) / ρ
    P = pressure(eos, ρ, e)
    return v, P
end

# ── Includes ───────────────────────────────────────────────────────────────────

include("../make_grid.jl")
include("../../fluxes/fluxes_include.jl")
include("../../solvers/HYDRO/FVM/Godunov/GodunovScheme.jl")
include("../../solvers/HYDRO/FVM/Godunov/GodunovSchemeDebug.jl")
include("../../limiters/minmod.jl")
include("../../limiters/vanleer.jl")
using Plots
using HDF5

# ── Simulation Parameters ──────────────────────────────────────────────────────

T = Float64
N = 1
origin = (0.0,)
domain = (2.0,)
nx = (720,)
ghost_zones = 3
cfl = 0.3
t_final = 2.0
h5_filename = "Test.h5"

states = [
    (density=1.0, velocity=(0.0,), pressure=1.0),
    (density=0.125, velocity=(0.0,), pressure=0.1)
]
walls = [1.0]

eos = IdealGasEOS(5/3)

# ── Grid ───────────────────────────────────────────────────────────────────────

grid_u = make_grid(:uniform, origin, domain, Int32.(nx), Int32(ghost_zones), mode=Verbose())

println("Uniform Grid:")
println("Centers: ", grid_u.axes[1].centers)
println("Faces:   ", grid_u.axes[1].faces)
println("dx:      ", grid_u.axes[1].dx)
println("Bounds:  ", grid_u.axes[1].bounds)
println("Active zones: ", grid_u.active_zones[1])
println()

# ── Initial Conditions ─────────────────────────────────────────────────────────

function ConstructShockTubeIC(grid, wall_positions::Vector{T}, states::Vector,
                              eos::AbstractEOS; mode::Verbosity = Standard()) where {T}

    log(mode, "Constructing Shock Tube initial conditions...")

    N = length(grid.axes)
    total = grid.total_zones[1]

    @assert length(states) == length(wall_positions) + 1

    walls = sort(wall_positions)

    log(mode, "Initializing variable centers...")
    density_c  = zeros(T, total)
    pressure_c = zeros(T, total)
    velocity_c = ntuple(_ -> zeros(T, total), N)

    log(mode, "Initializing variable faces...")
    density_f  = zeros(T, total)
    pressure_f = zeros(T, total)
    velocity_f = ntuple(_ -> zeros(T, total), N)

    xcenters = grid.axes[1].centers

    @inline function region_index(x)
        @inbounds for i in eachindex(walls)
            if x < walls[i]
                return i
            end
        end
        return length(walls) + 1
    end

    log(mode, "Assigning initial conditions to centers...")

    @inbounds for idx in 1:total
        r = region_index(xcenters[idx])
        s = states[r]

        density_c[idx]  = s.density
        pressure_c[idx] = s.pressure

        for d in 1:N
            velocity_c[d][idx] = s.velocity[d]
        end
    end

    log(mode, "Creating Primitive Fields...")

    centers_struct = PrimitiveFields{T,N}(density_c, velocity_c, pressure_c)
    faces_struct   = PrimitiveFields{T,N}(density_f, velocity_f, pressure_f)

    log(mode, "Returning Primitive Variables struct...")

    return PrimitiveVariables{T,N}(centers_struct, faces_struct)
end

function ConstructShuOsherIC(grid, eos::AbstractEOS;
                              mode::Verbosity = Standard())

    log(mode, "Constructing Shu-Osher initial conditions...")

    N     = length(grid.axes)
    total = grid.total_zones[1]
    T     = Float64

    density_c  = zeros(T, total)
    pressure_c = zeros(T, total)
    velocity_c = ntuple(_ -> zeros(T, total), N)

    density_f  = zeros(T, total)
    pressure_f = zeros(T, total)
    velocity_f = ntuple(_ -> zeros(T, total), N)

    xcenters = grid.axes[1].centers

    rho_L = T(3.857143)
    u_L   = T(2.629369)
    p_L   = T(10.33333)

    @inbounds for idx in 1:total
        x = xcenters[idx]

        if x < T(1.0)
            density_c[idx]     = rho_L
            pressure_c[idx]    = p_L
            velocity_c[1][idx] = u_L
        else
            density_c[idx]     = T(3.0) + T(0.2) * exp(x) * sin(T(100.0) * x)
            pressure_c[idx]    = T(1.0)
            velocity_c[1][idx] = T(0.0)
        end

        for d in 2:N
            velocity_c[d][idx] = T(0.0)
        end
    end

    centers_struct = PrimitiveFields{T,N}(density_c, velocity_c, pressure_c)
    faces_struct   = PrimitiveFields{T,N}(density_f, velocity_f, pressure_f)

    return PrimitiveVariables{T,N}(centers_struct, faces_struct)
end

# ── Construct Initial State ────────────────────────────────────────────────────

W = ConstructShockTubeIC(grid_u, walls, states, eos, mode = Verbose())
# W = ConstructShuOsherIC(grid_u, eos, mode = Standard())

E_centers = [total_energy(eos, W.centers.density[i], W.centers.velocity[1][i], W.centers.pressure[i])
             for i in eachindex(W.centers.density)]

E_faces   = [total_energy(eos, W.faces.density[i], W.faces.velocity[1][i], W.faces.pressure[i])
             for i in eachindex(W.faces.density)]

centers = ConservativeFields{Float64, 1}(
    W.centers.density,
    (W.centers.density .* W.centers.velocity[1],),
    E_centers
)

faces = ConservativeFields{Float64, 1}(
    W.faces.density,
    (W.faces.density .* W.faces.velocity[1],),
    E_faces
)

U = ConservativeVariables{Float64, 1}(centers, faces)
F = FluxVariables{Float64, 1}(
    zeros(length(grid_u.axes[1].faces)),
    (zeros(length(grid_u.axes[1].faces)),),
    zeros(length(grid_u.axes[1].faces))
)

include("../boundary_conditions/dispatch_boundary_condition.jl")
include("../boundary_conditions/reflecting.jl")

# ── Timestep ───────────────────────────────────────────────────────────────────

cs_init = [sound_speed(eos, W.centers.density[i], specific_energy(eos, W.centers.density[i], W.centers.pressure[i]))
           for i in eachindex(W.centers.density)]

dt = cfl * minimum(grid_u.axes[1].dx) / maximum(abs.(W.centers.velocity[1]) .+ cs_init)

println("Time step: ", dt, " s")

# ── Main Loop ──────────────────────────────────────────────────────────────────

t       = 0.0
counter = 0

println(length(U.centers.density))
plot(grid_u.axes[1].centers, U.centers.density, label="Density")

while t < t_final
    apply_reflecting_boundaries!(U, ghost_zones, nx[1])

    GodunovStep!(W, U, F, :Cubic, :vanleer, true, true, :Reflecting, :HLLC, eos, grid_u.axes[1].dx, dt, grid_u.total_zones[1], nx[1], ghost_zones, grid_u.axes[1].centers, :Standard)

    # Flux update (still needs GodunovStep! internals refactored — placeholder matches original)
    F.density_flux[1:end-1]      .= U.centers.momentum[1]
    F.momentum_flux[1][1:end-1]  .= U.centers.momentum[1].^2 ./ U.centers.density .+
                                     (eos.γ - 1) .* (U.centers.total_energy .-
                                     0.5 .* U.centers.momentum[1].^2 ./ U.centers.density)
    F.total_energy_flux[1:end-1] .= (U.centers.total_energy .+
                                     (eos.γ - 1) .* (U.centers.total_energy .-
                                     0.5 .* U.centers.momentum[1].^2 ./ U.centers.density)) .*
                                     (U.centers.momentum[1] ./ U.centers.density)

    # Primitive recovery — EOS-agnostic
    for i in eachindex(U.centers.density)
        W.centers.velocity[1][i], W.centers.pressure[i] = prim_from_cons(
            eos,
            U.centers.density[i],
            U.centers.momentum[1][i],
            U.centers.total_energy[i]
        )
    end
    W.centers.density .= U.centers.density

    # Floor
    for i in eachindex(W.centers.density)
        W.centers.density[i]  < 0.0 && (W.centers.density[i]  = 1e-10)
        W.centers.pressure[i] < 0.0 && (W.centers.pressure[i] = 1e-10)
    end

    # Timestep update
    cs = [sound_speed(eos, W.centers.density[i], specific_energy(eos, W.centers.density[i], W.centers.pressure[i]))
          for i in eachindex(W.centers.density)]

    global dt = cfl * grid_u.axes[1].dx / maximum(abs.(W.centers.velocity[1]) .+ cs)
    dt = min(dt, t_final - t)
    global t += dt

    # HDF5 output
    if counter % 10 == 0
        groupname = "step_$(counter)"
        println("Saving snapshot to $h5_filename in group $groupname")
        preset = isfile(h5_filename) ? "r+" : "w"
        h5open(h5_filename, preset) do file
            grp = create_group(file, groupname)
            grp["x"]             = grid_u.axes[1].centers
            grp["W/density"]     = W.centers.density
            grp["W/velocity"]    = W.centers.velocity[1]
            grp["W/pressure"]    = W.centers.pressure
            grp["U/momentum"]    = U.centers.momentum[1]
            grp["U/total_energy"] = U.centers.total_energy
            grp["time"]          = t
        end
    end

    global counter += 1
end

plot!(grid_u.axes[1].centers, U.centers.density, label="Density")