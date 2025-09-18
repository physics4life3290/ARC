




function Dispatch_Solver(input::LevelI;
    u::Union{Vector{Float64}, Nothing}=nothing
)
    ρ = input.state.ρ
    p = input.state.p
    γ = input.state.γ
    _grid = input.state.grid_points
    #println("Grid points are", _grid)
    coordinate_system = input.coordinate_system  
    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    cfl = 0.5
    solver = :GodunovScheme
    γ = 5/3
    mode = :Standard
    features = Symbol[]  # empty vector
    ghost_zones = 3
    zones = N - 2*ghost_zones
    boundary_condition = :None
    reconstruction = :Cubic
    limiter = :VanLeer
    flattening = true
    steepening = true
    riemanntype = :Exact
    operator_splitting = :Lie

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
    #GodunovStep!(W, U, F, reconstruction, limiter, flattening, steepening, boundary_condition, riemanntype, γ, spacing, dt, cfl, mode, features, N, zones, ghost_zones, _grid)
    
    HYDRO_Step!(W, U, F, dt,
        ghost_zones, N, spacing, zones, _grid,
        mode, features, cfl, solver,
        reconstruction, limiter, flattening, steepening,
        boundary_condition, riemanntype, γ
    )
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)

    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/N) .* F.density_flux
        mom_source = (1/N) .* F.momentum_flux
        tot_energy_source = (1/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/N) .* F.density_flux
        mom_source = (2/N) .* F.momentum_flux
        tot_energy_source = (2/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if operator_splitting == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        HYDRO_Step!(W, U, F, dt,
            ghost_zones, N, spacing, zones, _grid,
            mode, features, cfl, solver,
            reconstruction, limiter, flattening, steepening,
            boundary_condition, riemanntype, γ
        )
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Solver(input::LevelII;
    u::Union{Vector{Float64}, Nothing}=nothing
)
    ρ = input.state.ρ
    p = input.state.p
    γ = input.state.γ
    _grid = input.state.grid_points
    coordinate_system = input.coordinate_system
    N = length(ρ)
    spacing = _grid[2] - _grid[1]
    ghost_zones = 3
    mode = :Standard
    features = Symbol[]  # empty vector
    solver = input.solver
    cfl = input.Cfl
    boundary_condition = input.boundary_condition

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
    HYDRO_Step!(W, U, F, dt,
        ghost_zones, N, spacing, zones, _grid,
        mode, features, cfl, solver,
        reconstruction, limiter, flattening, steepening,
        boundary_condition, riemanntype, γ
    )
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)

    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/N) .* F.density_flux
        mom_source = (1/N) .* F.momentum_flux
        tot_energy_source = (1/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/N) .* F.density_flux
        mom_source = (2/N) .* F.momentum_flux
        tot_energy_source = (2/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if operator_splitting == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        HYDRO_Step!(W, U, F, dt,
            ghost_zones, N, spacing, zones, _grid,
            mode, features, cfl, solver,
            reconstruction, limiter, flattening, steepening,
            boundary_condition, riemanntype, γ
        )
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Solver(input::LevelIII;
    u::Union{Vector{Float64}, Nothing}=nothing
)
    ρ = input.state.ρ
    p = input.state.p
    γ = input.state.γ
    _grid = input.state.grid_points
    coordinate_system = input.coordinate_system
    boundary_condition = input.boundary_condition
    mode = input.mode
    features = input.features
    solver = input.solver
    cfl = input.Cfl
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
    HYDRO_Step!(W, U, F, dt,
        ghost_zones, N, spacing, zones, _grid,
        mode, features, cfl, solver,
        reconstruction, limiter, flattening, steepening,
        boundary_condition, riemanntype, γ
    )
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)

    if operator_splitting == :Strang
        dt = dt*2
    end

    if coordinate_system == :cylindrical
        dens_source = (1/N) .* F.density_flux
        mom_source = (1/N) .* F.momentum_flux
        tot_energy_source = (1/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    elseif coordinate_system == :spherical
        dens_source = (2/N) .* F.density_flux
        mom_source = (2/N) .* F.momentum_flux
        tot_energy_source = (2/N) .* F.total_energy_flux
        U.density_centers .-= dt/spacing .* dens_source
        U.momentum_centers .-= dt/spacing .* mom_source
        U.total_energy_centers .-= dt/spacing .* tot_energy_source
    end

    if operator_splitting == :Strang
        dt = dt/2
        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2
        F.density_flux .= U.momentum_centers
        F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
        F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)
        HYDRO_Step!(W, U, F, dt,
            ghost_zones, N, spacing, zones, _grid,
            mode, features, cfl, solver,
            reconstruction, limiter, flattening, steepening,
            boundary_condition, riemanntype, γ
        )
        dt = dt * 2
    end

    ρ .= U.density_centers
    u .= U.momentum_centers ./ U.density_centers
    p .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)

    return ρ, u, p
end

function Dispatch_Solver(input::LevelIV;
    u::Union{Vector{Float64}, Nothing}=nothing
)
    ρ = input.state.ρ
    p = input.state.p
    γ = input.state.γ
    _grid = input.state.grid_points
    coordinate_system = input.coordinate_system
    boundary_condition = input.boundary_condition
    mode = input.mode
    features = input.features
    solver = input.solver
    cfl = input.Cfl
    reconstruction = input.reconstruction
    limiter = input.limiter
    flattening = input.flattening
    steepening = input.steepening
    riemanntype = input.riemann_solver
    operator_splitting = :Strang
    N = length(ρ)
    spacing = _grid[2] - _grid[1]
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

    if operator_splitting == :Strang
        dt = dt/2
    end

    if operator_splitting == :None
        HYDRO_Step!(W, U, F, dt,
            ghost_zones, N, spacing, zones, _grid,
            mode, features, cfl, solver,
            reconstruction, limiter, flattening, steepening,
            boundary_condition, riemanntype, γ
        )
        return ρ, u, p
    end
    
    HYDRO_Step!(W, U, F, dt,
        ghost_zones, N, spacing, zones, _grid,
        mode, features, cfl, solver,
        reconstruction, limiter, flattening, steepening,
        boundary_condition, riemanntype, γ
    )
    
    F.density_flux .= U.momentum_centers
    F.momentum_flux .= U.momentum_centers.^2 ./ U.density_centers + (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)
    F.total_energy_flux .= (U.total_energy_centers .+ (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.momentum_centers .^ 2 ./ U.density_centers)) .* (U.momentum_centers ./ U.density_centers)

    if operator_splitting == :Strang
        dt = dt*2
    end

    if operator_splitting == :Lie || operator_splitting ==:Strang
        if coordinate_system == :cylindrical
            dens_source = (1/N) .* F.density_flux
            mom_source = (1/N) .* F.momentum_flux
            tot_energy_source = (1/N) .* F.total_energy_flux
            U.density_centers .-= dt/spacing .* dens_source
            U.momentum_centers .-= dt/spacing .* mom_source
            U.total_energy_centers .-= dt/spacing .* tot_energy_source
        elseif coordinate_system == :spherical
            dens_source = (2/N) .* F.density_flux
            mom_source = (2/N) .* F.momentum_flux
            tot_energy_source = (2/N) .* F.total_energy_flux
            U.density_centers .-= dt/spacing .* dens_source
            U.momentum_centers .-= dt/spacing .* mom_source
            U.total_energy_centers .-= dt/spacing .* tot_energy_source
        end
    end

    if operator_splitting == :Strang || operator_splitting ==:Lie
        ρ .= U.density_centers
        u .= U.momentum_centers ./ U.density_centers
        p .= (γ - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        return ρ, u, p
    end

end