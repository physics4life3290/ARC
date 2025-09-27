module ARC

using CUDA
using HDF5
using Plots
using ArgParse
using Dates
using Printf


#include("ARC_Include.jl")

#function CalculateFlux!(W, U, F)
#    F.density_flux .= U.momentum_centers 
#    F.momentum_flux .= U.momentum_centers .* W.velocity_centers .+ W.pressure_centers
#    F.total_energy_flux .= W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)
#end

include("../UI/initiate_UI.jl")
include("../UI/Tools/collect_primary_input.jl")
include("../UI/utility/prompt_setup.jl")
export run_Codex_Trials

struct PrimaryInput
    problem::Symbol           # :shocktube, :blastwave, etc.
    mode::Symbol              # :interactive, :batch
    features::Vector{Symbol}  # e.g., [:gravity, :cooling]
    dimension::Int            # 1, 2, or 3
    coordinate_system::Symbol # :cartesian, :cylindrical, :spherical
    solver::Symbol            # e.g., :godunov
    boundary_condition::Symbol # :reflecting, :outflow
    filename::String
end

PROBLEMS   = [:ShockTube, :BlastWave, :Custom]
SOLVERS    = [:FTCS, :LaxFriedrichs, :Richtmyer, :GodunovScheme]
DIMENSIONS = [1, 2, 3]  # using Symbols for type stability
COORDS     = [:Cartesian, :Spherical, :Cylindrical]
MODES      = [:Standard, :Parallel, :GPU, :HPC]
FEATURES   = [:Animate, :Benchmark, :Debug, :Verbose, :None]
BOUNDARIES = [:Reflecting, :Periodic, :Outflow, :None]

# --- Prompt functions ---
function prompt_choice(prompt::AbstractString, options::Vector{Symbol})
    while true
        println(prompt)
        for (i, opt) in enumerate(options)
            println("$(i). $opt")
        end
        print("> ")
        choice = tryparse(Int, readline())
        if choice !== nothing && 1 <= choice <= length(options)
            return options[choice]
        else
            println("Invalid choice. Enter a number between 1 and $(length(options)).")
        end
    end
end

function prompt_choice(prompt::AbstractString, options::Vector{Int})
    while true
        println(prompt)
        for (i, opt) in enumerate(options)
            println("$(i). $opt")
        end
        print("> ")
        choice = tryparse(Int, readline())
        if choice !== nothing && 1 <= choice <= length(options)
            return options[choice]
        else
            println("Invalid choice. Enter a number between 1 and $(length(options)).")
        end
    end
end

function prompt_multiple_choices(prompt::AbstractString, options::Vector{Symbol})
    while true
        println(prompt)
        for (i, opt) in enumerate(options)
            println("$(i). $opt")
        end
        println("Enter numbers separated by spaces (e.g., 1 3 4):")
        print("> ")
        input = readline()
        indices = try
            parse.(Int, split(input))
        catch
            println("Invalid input. Enter numbers separated by spaces.")
            continue
        end
        if all(1 .<= indices .<= length(options))
            return unique(options[indices])
        else
            println("Some choices are out of range. Try again.")
        end
    end
end

# --- Collect primary input ---
function collect_primary_input()
    println("=== Simulation Setup ===")

    problem_choice   = prompt_choice("Choose an example problem:", PROBLEMS)
    if problem_choice == :Custom
        @error "Custom problem setup is not yet implemented."
        return nothing
    end

    mode_choice      = prompt_choice("Choose a performance mode:", MODES)
    feature_choice   = prompt_multiple_choices("Choose additional features:", FEATURES)
    dimension_choice = prompt_choice("Choose simulation dimension:", DIMENSIONS)
    coord_choice     = prompt_choice("Choose coordinate system:", COORDS)
    solver_choice    = prompt_choice("Choose numerical solver:", SOLVERS)

    # Boundary conditions
    boundary_choice = (:Benchmark in feature_choice) ? :None :
                      prompt_choice("Choose boundary condition:", BOUNDARIES)

    println("Enter the output filename (without extension):")
    filename = readline()

    

    return PrimaryInput(
        problem_choice,
        mode_choice,
        feature_choice,
        dimension_choice,
        coord_choice,
        solver_choice,
        boundary_choice,
        filename
    )
end

# Prompt for a Float64 domain length
function prompt_domain_length(msg::String="Enter domain length:")
    println(msg)
    return parse(Float64, readline())
end

# Prompt for integer zones
function prompt_zones(msg::String="Enter number of zones:")
    println(msg)
    return parse(Int, readline())
end

# Prompt for ghost zones based on solver type
function prompt_ghost_zones(solver::Symbol)
    if solver in (:FTCS, :LaxFriedrichs, :Richtmyer)
        println("Enter number of ghost zones (≥ 1):")
        return parse(Int, readline())
    elseif solver == :GodunovScheme
        println("Enter number of ghost zones (1 for 1st order, 2 for MUSCL, 3 for PPM):")
        return parse(Int, readline())
    else
        error("Unsupported solver: $solver")
    end
end

struct GridInput{N}
    domain::NTuple{N,Float64}
    grid_center::NTuple{N,Float64}
    zones::NTuple{N,Int}
    ghost_zones::Int
    total_zones::NTuple{N,Int}
    coord_min::NTuple{N,Float64}
    coord_max::NTuple{N,Float64}
end

# --- Helper functions ---
function prompt_domain_length(msg::String="Enter domain length:")
    println(msg)
    return parse(Float64, readline())
end

function prompt_zones(msg::String="Enter number of zones:")
    println(msg)
    return parse(Int, readline())
end

function prompt_ghost_zones(solver::Symbol)
    if solver in (:FTCS, :LaxFriedrichs, :Richtmyer)
        println("Enter number of ghost zones (≥ 1):")
        return parse(Int, readline())
    elseif solver == :GodunovScheme
        println("Enter number of ghost zones (1 for 1st order, 2 for MUSCL, 3 for PPM):")
        return parse(Int, readline())
    else
        error("Unsupported solver: $solver")
    end
end

function parse_point(str::String)::Tuple{Vararg{Float64}}
    s = replace(str, r"[()]" => "")
    parts = split(strip(s), r"[\s,]+")
    parse.(Float64, parts) |> Tuple
end

# --- Main grid input collector ---
function collect_grid_input(primary_input::PrimaryInput)
    dim = primary_input.dimension
    cs  = primary_input.coordinate_system
    solver = primary_input.solver
    ghost_zones = nothing 
    
    if cs == :Cartesian
        # Cartesian grids
        if dim == 1
            domain = (prompt_domain_length("Enter domain length for 1D grid:"),)
            println("Enter grid center coordinate:")
            center = (parse(Float64, readline()),)
            zones = (prompt_zones(),)
        elseif dim == 2
            domain = (prompt_domain_length("Enter x-domain length:"), prompt_domain_length("Enter y-domain length:"))
            println("Enter grid center (x y) or (x,y):")
            center = parse_point(readline())
            zones = (prompt_zones("Number of x zones:"), prompt_zones("Number of y zones:"))
        elseif dim == 3
            domain = (
                prompt_domain_length("Enter x-domain length:"),
                prompt_domain_length("Enter y-domain length:"),
                prompt_domain_length("Enter z-domain length:")
            )
            println("Enter grid center (x y z) or (x,y,z):")
            center = parse_point(readline())
            zones = (
                prompt_zones("Number of x zones:"),
                prompt_zones("Number of y zones:"),
                prompt_zones("Number of z zones:")
            )
        else
            error("Unsupported dimension: $dim")
        end

        total_zones = ntuple(i -> zones[i] + 2*ghost_zones, dim)
        coord_min = ntuple(i -> -domain[i]/2 + center[i], dim)
        coord_max = ntuple(i ->  domain[i]/2 + center[i], dim)

        println("Total zones including ghost zones: ", prod(total_zones))

        return GridInput(domain, center, zones, ghost_zones, total_zones, coord_min, coord_max)

    elseif cs in (:Spherical, :Cylindrical)
        # Spherical/Cylindrical grids
        if dim == 1
            domain = (prompt_domain_length("Enter radial domain length:"),)
            println("Enter grid center (minimum coordinate):")
            center = (parse(Float64, readline()),)
            zones = (prompt_zones(),)
            total_zones = (zones[1] + 2*ghost_zones,)
            coord_min = (1e-12 + center[1],)
            coord_max = (domain[1] + center[1],)
        elseif dim == 2
            println("2D Spherical/Cylindrical not yet implemented")
            return nothing
        elseif dim == 3
            println("3D Spherical/Cylindrical not yet implemented")
            return nothing
        else
            error("Unsupported dimension: $dim")
        end

        ghost_zones = prompt_ghost_zones(solver)
        
        println("Total zones including ghost zones: ", prod(total_zones))
        return GridInput(domain, center, zones, ghost_zones, total_zones, coord_min, coord_max)
    else
        error("Unsupported coordinate system: $cs")
    end
end

function run_Codex_Trials()

    primary_input = collect_primary_input()
    grid_input = collect_grid_input(primary_input)
    println(primary_input)
    println(grid_input)
end

end