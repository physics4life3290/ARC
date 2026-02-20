




include("../make_grid.jl")

include("../../fluxes/fluxes_include.jl")
include("../../solvers/HYDRO/FVM/Godunov/GodunovScheme.jl")
#include("../boundary_conditions/reflecting.jl")
include("../../limiters/minmod.jl")
include("../../limiters/vanleer.jl")
using Plots
using HDF5
# ----------------------
# Example usage
# ----------------------
T = Float64          # numerical precision
N = 1                # number of dimensions (integer)
origin = (0.0,)      # physical coordinates (floats)
domain = (1.0,)      # physical size
nx = (200,)          # number of zones along each axis (integer)
ghost_zones = 3      # number of ghost zones (integer)
γ = 5/3              # physical constant (float)
cfl = 0.3            # CFL number (float)
t_final = 0.5        # simulation time (float)
h5_filename = "test.h5"

grid_u = make_grid(:uniform, origin, domain, Int32.(nx), Int32(ghost_zones))

println("Uniform Grid:")
println("Centers: ", grid_u.axes[1].centers)
#println("Centers 2: ", grid_u.axes[2].centers)
println("Faces:   ", grid_u.axes[1].faces)
println("dx:      ", grid_u.axes[1].dx)
println("Bounds:  ", grid_u.axes[1].bounds)
println("Active zones: ", grid_u.active_zones[1])
println()

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

function ConstructShockTubeIC(grid, wall_positions::Vector{T}, states::Vector, gamma::T) where {T}

    N = length(grid.axes)
    total = grid.total_zones[1]   # <--- cast to Int

    @assert length(states) == length(wall_positions) + 1

    walls = sort(wall_positions)

    density_c = zeros(T, total)
    pressure_c = zeros(T, total)
    velocity_c = ntuple(_ -> zeros(T, total), N)

    density_f = zeros(T, total)
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

    @inbounds for idx in 1:total
        r = region_index(xcenters[idx])
        s = states[r]

        density_c[idx] = s.density
        pressure_c[idx] = s.pressure

        for d in 1:N
            velocity_c[d][idx] = s.velocity[d]
        end
    end

    centers_struct = PrimitiveFields{T,N}(density_c, velocity_c, pressure_c)
    faces_struct   = PrimitiveFields{T,N}(density_f, velocity_f, pressure_f)

    return PrimitiveVariables{T,N}(centers_struct, faces_struct)
end

states = [
    (density=1.0, velocity=(0.0,), pressure=1.0),
    (density=0.125, velocity=(0.0,), pressure=0.1)
]

W = ConstructShockTubeIC(grid_u, [0.5], states, 1.4)

centers = ConservativeFields{Float64, 1}(W.centers.density, (W.centers.density .* W.centers.velocity[1], ), W.centers.pressure ./ (γ - 1) + 0.5 * (W.centers.density .* W.centers.velocity[1].^2))
faces = ConservativeFields{Float64, 1}(W.faces.density, (W.faces.density .* W.faces.velocity[1], ), W.faces.pressure ./ (γ - 1) + 0.5 * (W.faces.density .* W.faces.velocity[1].^2))

U = ConservativeVariables{Float64, 1}(centers, faces)
F = FluxVariables{Float64, 1}(zeros(length(grid_u.axes[1].faces)), (zeros(length(grid_u.axes[1].faces)), ), zeros(length(grid_u.axes[1].faces)))

include("../boundary_conditions/dispatch_boundary_condition.jl")
include("../boundary_conditions/reflecting.jl")

c = sqrt(γ * maximum(W.centers.pressure) / maximum(W.centers.density))
dt = cfl * minimum(grid_u.axes[1].dx) / (maximum(abs.(W.centers.velocity[1])) + c)

println("Time step: ", dt, " s")

t = 0.0
counter = 0
println(length(U.centers.density))
plot(grid_u.axes[1].centers, U.centers.density, label="Density")
while t < t_final 
    apply_reflecting_boundaries!(U, ghost_zones, nx[1])
    
    
    GodunovStep!(W, U, F, :Cubic, :vanleer, true, true, :Reflecting, :Exact, γ, grid_u.axes[1].dx, dt, grid_u.total_zones[1], nx[1], ghost_zones, grid_u.axes[1].centers, :Standard)

    F.density_flux[1:end-1] .= U.centers.momentum[1]
    F.momentum_flux[1][1:end-1] .= U.centers.momentum[1].^2 ./ U.centers.density + (γ - 1) .* (U.centers.total_energy .- 0.5 .* U.centers.momentum[1] .^ 2 ./ U.centers.density)
    F.total_energy_flux[1:end-1] .= (U.centers.total_energy .+ (γ - 1) .* (U.centers.total_energy .- 0.5 .* U.centers.momentum[1] .^ 2 ./ U.centers.density)) .* (U.centers.momentum[1] ./ U.centers.density)

    W.centers.density .= U.centers.density
    W.centers.velocity[1] .= U.centers.momentum[1] ./ U.centers.density
    W.centers.pressure .= (γ - 1) .* (U.centers.total_energy .- 0.5 .* U.centers.density .* W.centers.velocity[1] .^ 2)

    
    for i in 1:length(W.centers.density)

        if W.centers.density[i] < 0.0
            W.centers.density[i] = 1e-10
        end
        if W.centers.pressure[i] < 0.0
            W.centers.pressure[i] = 1e-10
        end

    end
    
    global c = sqrt.(γ .* W.centers.pressure ./ W.centers.density)
    global dt = cfl * grid_u.axes[1].dx / maximum(abs.(W.centers.velocity[1]) .+ c)
    dt = min(dt, t_final - t)
    global t += dt

    if counter % 10 == 0

        global groupname = "step_$(counter)"
        println("Saving snapshot to $h5_filename in group $groupname")
        preset = isfile(h5_filename) ? "r+" : "w"
        h5open(h5_filename, preset) do file  # "a" = append mode
            grp = create_group(file, groupname)

            grp["x"] = grid_u.axes[1].centers

            grp["W/density"]      = W.centers.density
            grp["W/velocity"]     = W.centers.velocity[1]
            grp["W/pressure"]     = W.centers.pressure

            grp["U/momentum"]     = U.centers.momentum[1]
            grp["U/total_energy"] = U.centers.total_energy

            grp["time"] = t

        end

    end

    global counter += 1
end
plot!(grid_u.axes[1].centers, U.centers.density, label="Density")