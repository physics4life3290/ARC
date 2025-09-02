




function Dispatch_LaxFriedrichs_I(
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
    LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, N, spacing, mode, features, cfl)

    return W, U
end


function Dispatch_LaxFriedrichs_II(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64}, 
    cfl::Float64, 
    γ::Float64;
    u::Union{Vector{Float64}, Nothing}=nothing
)

    ####### ENFORCEMENTS #######
    N = length(ρ)
    spacing = grid[2] - grid[1]
    mode = :Standard
    features = Symbol[]  # empty vector
    ghost_zones = 3

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)
    # Build primitive variable object
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u .^ 2, nothing, nothing, nothing)
    F = FluxVariables(zeros(length(ρ)), zeros(length(ρ)), zeros(length(ρ)))
    ####### CALL SOLVER #######
    LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, N, spacing, mode, features, cfl)
    return W, U
end

function Dispatch_LaxFriedrichs_III(
    ρ::Vector{Float64},
    p::Vector{Float64},
    grid::Vector{Float64}, 
    cfl::Float64, 
    γ::Float64, 
    ghost_zones::Int64,
    mode::Symbol, 
    features::Vector{Symbol};
    u::Union{Vector{Float64}, Nothing}=nothing
) 

    ####### ENFORCEMENTS #######
    N = length(ρ)
    spacing = grid[2] - grid[1]

    if u === nothing
        u = zeros(N)  # default velocity
    else
        @assert length(u) == N "Velocity array must match density array length"
    end

    # Compute sound speed and timestep
    c = sqrt.(γ .* p ./ ρ)
    dt = cfl * spacing / maximum(abs.(u) .+ c)

    ϵ = p ./ ((γ-1) .* ρ)
    # Build primitive variable object
    W = PrimitiveVariables(ρ, u, p, ϵ, nothing, nothing, nothing, nothing)
    U = ConservativeVariables(ρ, ρ .* u, ρ .* ϵ .+ 0.5 .* ρ .* u .^ 2, nothing, nothing, nothing)
    F = FluxVariables(zeros(length(ρ)), zeros(length(ρ)), zeros(length(ρ)))
    ####### CALL SOLVER #######
    LaxFriedrichs_Step!(W, U, F, dt, ghost_zones, N, spacing, mode, features, cfl)

    return W, U
end