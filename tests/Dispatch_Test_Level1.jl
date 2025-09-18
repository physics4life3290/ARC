include("../src/ARC_Include.jl")

using Plots

# -------------------------------
# Define the grid and base state. This is the Users input!
# -------------------------------
N = 1000
x = range(0.0000000001, stop=1.0000000001, length=N) |> collect
ρ0 = ones(N)
p0 = ones(N)
γ0 = 5/3
coordinate_system = :spherical
grid_center = 0.5000000001
domain = 1.0

# Initial discontinuity (shock tube setup)
for i in 1:N
    if x[i] < grid_center
        ρ0[i] = 1.0
        p0[i] = 1.0
    else
        ρ0[i] = 0.125
        p0[i] = 0.1
    end
end

u0 = zeros(N)  # stationary initial velocity

plot(x, ρ0, label="Density")
plot!(x, u0, label="Velocity")
plot!(x, p0, label="Pressure")
title!("1D Shock Tube Initial Conditions")

spacing = x[2] - x[1]
ghost_zones = 3 
boundary_condition = :Reflecting
operator_splitting = :Strang

# --------------------------
# Create temp grid for step 
#---------------------------
dispatch_grid = nothing
if coordinate_system == :Cartesian
    dispatch_grid = Construct1DCartesian(domain, N, ghost_zones, grid_center, "cm")
elseif coordinate_system == :cylindrical || coordinate_system == :spherical
    dispatch_grid = Construct1DSpherical(domain, N, ghost_zones, grid_center, "cm")
end
println(dispatch_grid)

# Extend intial conditions to include ghost zones
ρ0 = vcat(fill(ρ0[1], ghost_zones), ρ0, fill(ρ0[end], ghost_zones))
u0 = vcat(fill(u0[1], ghost_zones), u0, fill(u0[end], ghost_zones))
p0 = vcat(fill(p0[1], ghost_zones), p0, fill(p0[end], ghost_zones))

plot(dispatch_grid.coord1.all_centers, ρ0, label="Density")
plot!(dispatch_grid.coord1.all_centers, u0, label="Velocity")
plot!(dispatch_grid.coord1.all_centers, p0, label="Pressure")
title!("1D Shock Tube Initial Conditions")

# -------------------------------
# Time evolution parameters
# -------------------------------

cfl = 0.3
t = 0.0
t_final = 0.2
c = sqrt.(γ0 .* p0 ./ ρ0)
dt = cfl * spacing / maximum(abs.(u0) .+ c)
dt = min(dt, t_final - t)
t += dt
ghost_zones = 3

println("The maximum sound speed is: ", maximum(c))
println("The timestep is: ", dt)

# Initial solution fields
ρ = copy(ρ0)
u = copy(u0)
p = copy(p0)
ϵ = p ./ ((γ0-1) .* ρ)
# Construct Primitive, Conservative, and Flux variables
W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
U = ConservativeVariables(ρ, ρ .* u, ρ .* (ϵ .+ 0.5 .* u.^2) .+ p, nothing, nothing, nothing)
F = FluxVariables(ρ .* u, ρ .* u .^ 2 .+ p, u .* (ρ .* (ϵ .+ 0.5 .* u.^2) .+ p))


apply_boundary_conditions(boundary_condition, U, N, ghost_zones)

if operator_splitting == :Strang
    dt /= 2
end

HYDRO_Step!(W, U, F, dt, ghost_zones, N, spacing, N, x, :Standard, Symbol[], cfl, :GodunovScheme, :Cubic, :VanLeer, true, true, boundary_condition, :Exact, γ0)


if coordinate_system != :cartesian
    
    if operator_splitting == :Strang
        dt *= 2
    end

    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers .^ 2 ./ U.density_centers .+ (γ0 - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (γ0 - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
    #apply_boundary_conditions(boundary_condition, U, N, ghost_zones)
    if coordinate_system == :cylindrical
        dens_source = (1 ./ dispatch_grid.coord1.all_centers) .* F.density_flux
        mom_source = (1 ./ dispatch_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (1 ./ dispatch_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2 ./ dispatch_grid.coord1.all_centers) .* F.density_flux
        mom_source = (2 ./ dispatch_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (2 ./ dispatch_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if operator_splitting == :Strang
        dt /= 2
        apply_boundary_conditions(boundary_condition, U, N, ghost_zones)
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers .^ 2 ./ U.density_centers .+ (γ0 - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (γ0 - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        HYDRO_Step!(W, U, F, dt, ghost_zones, N, spacing, N, x, :Standard, Symbol[], cfl, :GodunovScheme, :Cubic, :VanLeer, true, true, boundary_condition, :Exact, γ0)
    end

end


plot(dispatch_grid.coord1.all_centers, U.density_centers, label="Density")
plot!(dispatch_grid.coord1.all_centers, U.momentum_centers, label="Momentum")
plot!(dispatch_grid.coord1.all_centers, U.total_energy_centers, label="Total Energy")
# This works so far!
