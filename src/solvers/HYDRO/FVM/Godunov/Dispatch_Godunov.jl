









function Dispatch_Godunov_I(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64};
    u::Union{Vector{Float64}, Nothing}=nothing
)
    N = length(ρ)
    spacing = grid[2] - grid[1]
    γ = 5/3
    mode = :Standard
    features = Symbol[]  # empty vector
    cfl = 0.5
    ghost_zones = 3
    zones = N - 2*ghost_zones
    boundary_condition = :None
    reconstruction = :Parabolic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact

    @assert length(p) == N "Pressure array must match density array length"

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)

    # Build primitive and conservative variable objects
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u.^2, nothing, nothing, nothing)
    F = FluxVariables(zeros(N), zeros(N), zeros(N))

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, )
    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end


function Dispatch_Godunov_II(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64}, 
    cfl::Float64, 
    γ::Float64;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = grid[2] - grid[1]
    mode = :Standard
    features = Symbol[]  # empty vector
    ghost_zones = 3
    zones = N - 2*ghost_zones
    boundary_condition = :None
    reconstruction = :Parabolic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact

    @assert length(p) == N "Pressure array must match density array length"

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)

    # Build primitive and conservative variable objects
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u.^2, nothing, nothing, nothing)
    F = FluxVariables(zeros(N), zeros(N), zeros(N))

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones)
    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Godunov_III(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64}, 
    cfl::Float64, 
    γ::Float64, 
    ghost_zones::Int64,
    mode::Symbol, 
    features::Vector{Symbol}, 
    boundary_condition::Symbol;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = grid[2] - grid[1]
    zones = N - 2*ghost_zones
    reconstruction = :Parabolic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact

    @assert length(p) == N "Pressure array must match density array length"

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)

    # Build primitive and conservative variable objects
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u.^2, nothing, nothing, nothing)
    F = FluxVariables(zeros(N), zeros(N), zeros(N))

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones)

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Godunov_IV(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64}, 
    cfl::Float64, 
    γ::Float64, 
    ghost_zones::Int64,
    mode::Symbol, 
    features::Vector{Symbol}, 
    boundary_condition::Symbol, 
    reconstruction, 
    limiter,
    flattening, 
    steepening,
    riemanntype;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = grid[2] - grid[1]
    zones = N - 2*ghost_zones

    @assert length(p) == N "Pressure array must match density array length"

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)

    # Build primitive and conservative variable objects
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u.^2, nothing, nothing, nothing)
    F = FluxVariables(zeros(N), zeros(N), zeros(N))

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones)

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end