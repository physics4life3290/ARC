




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

include("../make_grid.jl")

include("../../fluxes/fluxes_include.jl")
include("../../solvers/HYDRO/FVM/Godunov/GodunovScheme_gravity.jl")
include("../../solvers/HYDRO/FVM/Godunov/GodunovSchemeDebug.jl")
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
domain = (6.957e10,)      # physical size
nx = (1000,)          # number of zones along each axis (integer)
ghost_zones = 3      # number of ghost zones (integer)
γ = 5/3              # physical constant (float)
cfl = 0.3            # CFL number (float)
t_final = 1E8        # simulation time (float)
h5_filename = "test.h5"
G = 6.6743e-8

grid_u = make_grid(:uniform, origin, domain, Int32.(nx), Int32(ghost_zones))

println("Uniform Grid:")
println("Centers: ", grid_u.axes[1].centers)
#println("Centers 2: ", grid_u.axes[2].centers)
println("Faces:   ", grid_u.axes[1].faces)
println("dx:      ", grid_u.axes[1].dx)
println("Bounds:  ", grid_u.axes[1].bounds)
println("Active zones: ", grid_u.active_zones[1])
println()

# --- Fix 1: Matching N in PrimitiveFields ---
# The N in NTuple{N, Vector{T}} must match the N in PrimitiveFields{T,N}
function ConstructPolytropeIC(grid, ρc::T, γ::T) where T
    r = grid.axes[1].centers
    # Get N from the grid type or dimensionality
    # If this is 1D, N must be 1.
    D = 1 

    n = 1 / (γ - 1)
    R_star = maximum(r) * 0.9 
    
    density_c = [ρc * max(1 - (ri/R_star)^2, 0.0)^n for ri in r]
    pressure_c = [ρ^γ for ρ in density_c]
    
    # Velocity must be an NTuple of length D
    velocity_c = ntuple(_ -> zeros(T, length(r)), D)

    centers_struct = PrimitiveFields{T, D}(density_c, velocity_c, pressure_c)
    
    # Ensure faces_struct matches the same N
    faces_v = ntuple(_ -> zeros(T, length(grid.axes[1].faces)), D)
    faces_struct = PrimitiveFields{T, D}(
        zeros(T, length(grid.axes[1].faces)),
        faces_v,
        zeros(T, length(grid.axes[1].faces))
    )

    return PrimitiveVariables{T, D}(centers_struct, faces_struct)
end

# --- Fix 2: Enclosed Mass Logic ---
function compute_enclosed_mass(grid, density::Vector{T}) where T
    r_faces = grid.axes[1].faces
    n_cells = length(density)
    mass = zeros(T, n_cells)
    
    running_mass = 0.0
    for i in 1:n_cells
        r_inner = r_faces[i]
        r_outer = r_faces[i+1]
        shell_volume = 4/3 * π * (r_outer^3 - r_inner^3)
        running_mass += density[i] * shell_volume
        mass[i] = running_mass
    end
    return mass
end

# --- Fix 3: Gravitational Potential Array Handling ---
# Your original code had an indexing mismatch between active zones and the full grid
W = ConstructPolytropeIC(grid_u, 5.7, γ)
enclosed_mass = compute_enclosed_mass(grid_u, W.centers.density)

# Initialize potential for the full grid (including ghost zones)
gravity_full = zeros(T, length(grid_u.axes[1].centers))

# Compute potential using the full grid centers to avoid slicing errors
function compute_gravitational_potential(grid, enclosed_mass::Vector{T}, G::T) where T
    r = grid.axes[1].centers
    N = length(r)
    Φ = zeros(T, N)

    # Standard radial integration: dΦ/dr = GM/r^2
    for i in 2:N
        dr = r[i] - r[i-1]
        # Use average mass or midpoint for better accuracy
        Φ[i] = Φ[i-1] + (G * enclosed_mass[i] / r[i]^2) * dr
    end
    return Φ
end

gravity_full = compute_gravitational_potential(grid_u, enclosed_mass, G)
gravity_full[1:ghost_zones] = gravity_full[ghost_zones+1:2*ghost_zones]

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
    # 1. Update Primitive Variables from Conserved (Sync step)
    # This ensures Reconstruction in GodunovStep! uses the latest data
    @inbounds for i in 1:grid_u.total_zones[1]
        # Apply physical floors to prevent NaN in sqrt(p/ρ)
        U.centers.density[i] = max(U.centers.density[i], 1e-10)
        
        W.centers.density[i] = U.centers.density[i]
        W.centers.velocity[1][i] = U.centers.momentum[1][i] / U.centers.density[i]
        
        # Internal Energy = Total - Kinetic
        e_int = U.centers.total_energy[i] - 0.5 * U.centers.density[i] * W.centers.velocity[1][i]^2
        W.centers.pressure[i] = max((γ - 1.0) * e_int, 1e-10)
    end

    # 2. Update Gravity (Optional but recommended for dynamic polytropes)
    # enclosed_mass = compute_enclosed_mass(grid_u, W.centers.density)
    # gravity_full  = compute_gravitational_potential(grid_u, enclosed_mass, G)

    # 3. Apply Boundary Conditions
    apply_reflecting_boundaries!(U, ghost_zones, nx[1])
    
    # 4. The Solver Step (This now handles Fluxes + Lie Splitting Source Terms)
    GodunovStep!(
        W, U, F, gravity_full, 
        :Cubic, :vanleer, true, true, 
        :Reflecting, :HLLC, γ, 
        grid_u.axes[1].dx, dt, 
        grid_u.total_zones[1], nx[1], ghost_zones, 
        grid_u.axes[1].centers, :Standard
    )    

    # 5. Adaptive Timestep (CFL)
    # Recalculate sound speed 'c' based on the updated W
    c_soundspeed = sqrt.(γ .* W.centers.pressure ./ W.centers.density)
    max_wave_speed = maximum(abs.(W.centers.velocity[1]) .+ c_soundspeed)
    
    global dt = cfl * grid_u.axes[1].dx / max_wave_speed
    
    # Final step check
    if t + dt > t_final
        dt = t_final - t
    end

    # 6. Data Logging / HDF5
    if counter % 10 == 0
        global groupname = "step_$(counter)"
        println("t = $(round(t, digits=4)) | dt = $(round(dt, digits=6)) | Snapshot: $groupname")
        
        h5open(h5_filename, isfile(h5_filename) ? "r+" : "w") do file
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

    global t += dt
    global counter += 1
end

plot!(grid_u.axes[1].centers, U.centers.density, label="Density")
