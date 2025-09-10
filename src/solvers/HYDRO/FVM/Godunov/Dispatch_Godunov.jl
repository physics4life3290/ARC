









function Dispatch_Godunov_I(
    ρ::Vector{Float64},
    p::Vector{Float64},
    _grid::Vector{Float64}, 
    coordinate_system::Symbol;
    u::Union{Vector{Float64}, Nothing}=nothing
)
    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    cfl = 0.5
    γ = 5/3
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
    operator_splitting = :Strang

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

    if operator_splitting == :Strang
        dt = dt/2
    end

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
    
    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (1/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (1/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (2/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (2/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if user_input.Solver_Input.split_choice == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end


function Dispatch_Godunov_II(
    ρ::Vector{Float64},
    p::Vector{Float64},
    _grid::Vector{Float64},
    coordinate_system::Symbol, 
    cfl::Float64, 
    γ::Float64;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    ghost_zones = 3
    mode = :Standard
    features = Symbol[]  # empty vector
    boundary_condition = :None

    zones = N - 2*ghost_zones
    reconstruction = :Parabolic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact
    operator_splitting = :Strang

    @assert length(p) == N "Pressure arra  y must match density array length"

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

    if operator_splitting == :Strang
        dt = dt/2
    end

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
    
    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (1/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (1/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (2/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (2/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if user_input.Solver_Input.split_choice == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Godunov_III(
    ρ::Vector{Float64},
    p::Vector{Float64},
    _grid::Vector{Float64},
    coordinate_system::Symbol, 
    cfl::Float64, 
    γ::Float64, 
    ghost_zones::Int64,
    mode::Symbol, 
    features::Vector{Symbol}, 
    boundary_condition::Symbol;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    zones = N - 2*ghost_zones
    reconstruction = :Parabolic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact
    operator_splitting = :Strang

    @assert length(p) == N "Pressure arra  y must match density array length"

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

    if operator_splitting == :Strang
        dt = dt/2
    end

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
    
    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (1/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (1/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (2/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (2/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if user_input.Solver_Input.split_choice == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Godunov_IV(
    ρ::Vector{Float64},
    p::Vector{Float64},
    _grid::Vector{Float64}, 
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
    riemanntype, 
    operator_splitting;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    zones = N - 2*ghost_zones

    @assert length(p) == N "Pressure arra  y must match density array length"

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

    if operator_splitting == :Strang
        dt = dt/2
    end

    # Call solver
    GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
    
    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (1/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (1/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/_grid.coord1.all_centers) .* F.density_flux
        mom_source = (2/_grid.coord1.all_centers) .* F.momentum_flux
        tot_energy_source = (2/_grid.coord1.all_centers) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if user_input.Solver_Input.split_choice == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
    
    return ρ, u, p
end