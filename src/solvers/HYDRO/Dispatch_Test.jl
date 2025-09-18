


include("DispatchStructs.jl")
include("../../ARC_Include.jl")
gr()
# ----------------------------------------------------------------------
# Mock Input
# ----------------------------------------------------------------------
dens  = [1.0, 1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125, 0.125, 0.125]
press = [1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1]
gamma = 5/3
grid  = collect(0.0:0.1:1.0)

coord_sys = :Cartesian
BSTATE    = BaseState(dens, press, gamma, grid, coord_sys)
L1        = LevelI(BSTATE)

# ----------------------------------------------------------------------
# Solver Parameters
# ----------------------------------------------------------------------
cfl                = 0.3
t                  = 0.0
t_final            = 1.0
ghost_zones        = 3
zones              = length(dens)
total_zones        = zones + 2 * ghost_zones
spacing            = grid[2] - grid[1]

mode               = :Standard
features           = Symbol[]
solver             = :GodunovScheme
reconstruction     = :Parabolic
limiter            = :Vanleer
flat, steep        = true, true
riemanntype        = :Exact
boundary_condition = :None

# ----------------------------------------------------------------------
# Dispatch Grid and Variables
# ----------------------------------------------------------------------
dispatch_grid  = collect(-grid[1] - spacing*ghost_zones : spacing : grid[end] + spacing*ghost_zones)
dispatch_dens  = vcat(fill(dens[1],  ghost_zones), dens, fill(dens[end],  ghost_zones))
dispatch_press = vcat(fill(press[1], ghost_zones), press, fill(press[end], ghost_zones))

W = PrimitiveVariables(
        dispatch_dens,
        dispatch_press,
        zeros(length(dispatch_dens)),  # velocities
        zeros(length(dispatch_dens)), nothing, nothing, nothing, nothing
    )

U = ConservativeVariables(
        W.density_centers,
        W.density_centers .* W.velocity_centers,
        W.internal_energy_centers .+ 0.5 .* W.density_centers .* W.velocity_centers,
        nothing, nothing, nothing
    )

F = FluxVariables(
        U.momentum_centers,
        U.momentum_centers .* W.velocity_centers .+ W.pressure_centers,
        W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)
    )

# ----------------------------------------------------------------------
# Time Step
# ----------------------------------------------------------------------
for i in 1:2
    c  = sqrt.(gamma .* W.pressure_centers ./ W.density_centers)
    dt = cfl * spacing / maximum(abs.(W.velocity_centers) .+ c)
    dt = min(dt, t_final - t)   # assumes t_final and t are defined

    # ----------------------------------------------------------------------
    # Solver Step
    # ----------------------------------------------------------------------
    HYDRO_Step!(
        W, U, F, dt,
        ghost_zones, total_zones, spacing, zones, grid,
        mode, features, cfl, solver,
        reconstruction, limiter, flat, steep,
        boundary_condition, riemanntype, gamma
    )
end
plot(dispatch_grid[1:end-1], W.density_centers, label="Density", xlabel="x", ylabel="Density", legend=:topright)
plot!(dispatch_grid[1:end-1], W.pressure_centers, label="Pressure", xlabel="x", ylabel="Pressure", legend=:topright)
plot!(dispatch_grid[1:end-1], W.velocity_centers, label="Velocity", xlabel="x", ylabel="Velocity", legend=:topright)