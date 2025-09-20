include("../src/ARC_Include.jl")

#using Plots

# -------------------------------
# Define the grid and base state. This is the Users input!
# -------------------------------
N = 100
x = range(0.0000000000, stop=1.0000000000, length=N) |> collect
ρ0 = ones(N)
p0 = ones(N)
γ0 = 5/3
coordinate_system = :Cartesian
# Initial discontinuity (shock tube setup)
for i in 1:N
    if x[i] < (x[end]-x[1])/2
        ρ0[i] = 1.0
        p0[i] = 1.0
    else
        ρ0[i] = 0.125
        p0[i] = 0.1
    end
end

u0 = zeros(N)  # stationary initial velocity

println("Starting Level I Dispatch Test...")
input = LevelI(BaseState(ρ0, u0, p0, γ0, x, coordinate_system))


# This is where the Dispatch function begins
#input_extract_time1 = time_ns()
ρ, u, p, γ, x, coordinate_system = input.state.ρ, input.state.u, input.state.p, input.state.γ, input.state.grid_points, input.state.coordinate_system
#input_extract_time2 = time_ns()
#println("Input Extraction Time: ", (input_extract_time2 - input_extract_time1)/1e6, " ms")
dx = x[2] - x[1]
ghost_zones = 3
N = length(x)
boundary_condition = :Reflecting
operator_splitting = :Strang
cfl, t_final = 0.3, 0.2
domain = x[end] - x[1]
grid_center = 0.5 * (x[1] + x[end])

grid_creation_time1 = time_ns()
dispatch_grid = Construct1DCartesian(domain, N, ghost_zones, grid_center, "cm")
grid_creation_time2 = time_ns()
println("Grid Creation Time: ", (grid_creation_time2 - grid_creation_time1)/1e6, " ms")

function fill_ghosts!(dest::Vector{Float64}, src::Vector{Float64}, ghost_zones::Int)
    N = length(src)
    @inbounds begin
        for i in 1:ghost_zones
            dest[i] = src[1]
            dest[ghost_zones+N+i] = src[end]
        end
        for i in 1:N
            dest[ghost_zones+i] = src[i]
        end
    end
end

t0 = time_ns()
ng = N + 2*ghost_zones
ρ = Vector{Float64}(undef, ng)
u   = Vector{Float64}(undef, ng)
p   = Vector{Float64}(undef, ng)
ϵ   = Vector{Float64}(undef, ng)
fill_ghosts!(ρ, ρ0, ghost_zones)
fill_ghosts!(u, u0, ghost_zones)
fill_ghosts!(p, p0, ghost_zones)

@inbounds for i in 1:ng
    ϵ[i] = p[i] / ((γ0 - 1) * ρ[i])
end
t1 = time_ns()
println("Preallocation + Ghost Filling + ϵ Calculation: ", (t1-t0)/1e6, " ms")

W = PrimitiveVariables(ρ, u, p, ϵ, zeros(ng), zeros(ng), zeros(ng), zeros(ng))
U = ConservativeVariables(ρ, similar(ρ), similar(ρ), zeros(ng), zeros(ng), zeros(ng))
F = FluxVariables(similar(ρ), similar(ρ), similar(ρ))

t_struct_init0 = time_ns()
@inbounds for i in 1:ng
    U.momentum_centers[i] = ρ[i] * u[i]
    U.total_energy_centers[i] = ρ[i]*(ϵ[i] + 0.5*u[i]^2) + p[i]
    F.density_flux[i] = ρ[i] * u[i]
    F.momentum_flux[i] = ρ[i]*u[i]^2 + p[i]
    F.total_energy_flux[i] = u[i]*(ρ[i]*(ϵ[i] + 0.5*u[i]^2) + p[i])
end
t_struct_init1 = time_ns()
println("Struct Initialization Loops: ", (t_struct_init1 - t_struct_init0)/1e6, " ms")

function compute_timestep(u::Vector{Float64}, ρ::Vector{Float64}, p::Vector{Float64},
                          γ::Float64, dx::Float64, cfl::Float64, t_final::Float64)
    @inbounds begin
        max_speed = 0.0
        for i in eachindex(u)
            c = sqrt(γ * p[i] / ρ[i])
            speed = abs(u[i]) + c
            if speed > max_speed
                max_speed = speed
            end
        end
        dt = cfl * dx / max_speed
        return min(dt, t_final)
    end
end


t0 = time_ns()
dt = compute_timestep(u, ρ, p, γ0, dx, cfl, t_final)
t1 = time_ns()
println("dt Computation Time: ", (t1-t0)/1e6, " ms")

t0 = time_ns()
apply_boundary_conditions(boundary_condition, U, N, ghost_zones)
t1 = time_ns()
println("Boundary Condition Application Time: ", (t1-t0)/1e6, " ms")

# -------------------------------
# Optional Strang adjustment
# -------------------------------
if operator_splitting == :Strang
    dt /= 2
end

# -------------------------------
# Evolve solution
# -------------------------------
features = Vector{Symbol}(undef, 0)        # preallocate empty vector of known type
boundary_condition = :Reflecting           # ensure Symbol, not String
x = collect(Float64, x)                    # ensure Vector{Float64}
N = Int(N)
dx = Float64(dx)
γ0 = Float64(γ0)
cfl = Float64(cfl)
flatten, steepen = true, true
mode, solver, reconstruction, limiter, riemanntype = :Standard, :GodunovScheme, :Cubic, :VanLeer, :Exact

# Define a type for the solver selection
abstract type SolverMethod end
struct GodunovSolver <: SolverMethod end
struct RichtmyerSolver <: SolverMethod end
struct LaxFriedrichsSolver <: SolverMethod end

# Dispatcher function
function step!(solver::SolverMethod, W::PrimitiveVariables, U::ConservativeVariables,
               F::FluxVariables; kwargs...)
    # Call the appropriate solver based on type
    _step!(solver, W, U, F; kwargs...)
end

# Godunov
function _step!(::GodunovSolver, W, U, F; kwargs...)
    GodunovStep!(W, U, F,
        kwargs[:reconstruction], kwargs[:limiter], kwargs[:flattening], kwargs[:steepening],
        kwargs[:boundary_condition], kwargs[:riemanntype], kwargs[:γ],
        kwargs[:spacing], kwargs[:dt], kwargs[:cfl], kwargs[:mode], kwargs[:features],
        kwargs[:total_zones], kwargs[:zones], kwargs[:ghost_zones], kwargs[:grid_points])
end

# Richtmyer
function _step!(::RichtmyerSolver, W, U, F; kwargs...)
    RichtmyerStep!(W, U, F,
        kwargs[:ghost_zones], kwargs[:total_zones], kwargs[:spacing],
        kwargs[:boundary_condition], kwargs[:cfl], kwargs[:dt], kwargs[:γ],
        kwargs[:mode], kwargs[:features], kwargs[:zones])
end

# Lax-Friedrichs
function _step!(::LaxFriedrichsSolver, W, U, F; kwargs...)
    LaxFriedrichs_Step!(W, U, F,
        kwargs[:dt], kwargs[:ghost_zones], kwargs[:total_zones], kwargs[:spacing],
        kwargs[:mode], kwargs[:features], kwargs[:cfl])
end

solver = GodunovSolver()  # or RichtmyerSolver(), LaxFriedrichsSolver()
step!(solver, W, U, F;
      reconstruction=:Cubic, limiter=:minmod, flattening=true, steepening=true,
      boundary_condition=:Reflecting, riemanntype=:Exact, γ=1.4,
      spacing=0.01, dt=1e-4, cfl=0.8, mode=:Standard,
      features=[:None], total_zones=100, zones=98, ghost_zones=2,
      grid_points=100)
