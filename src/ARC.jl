module ARC

using CUDA
using HDF5
using Plots
using ArgParse # For command line argument parsing access

const MaybeVector{T} = Union{Nothing, Vector{T}}

struct PrimitiveVariables 
    density_centers::Vector{Float64}
    velocity_centers::Vector{Float64}
    pressure_centers::Vector{Float64}
    internal_energy_centers::Vector{Float64}
    density_faces::MaybeVector{Float64}
    velocity_faces::MaybeVector{Float64}
    pressure_faces::MaybeVector{Float64}
    internal_energy_faces::MaybeVector{Float64}
end

struct ConservativeVariables
    density_centers::Vector{Float64}
    momentum_centers::Vector{Float64}
    total_energy_centers::Vector{Float64}
    density_faces::MaybeVector{Float64}
    momentum_faces::MaybeVector{Float64}
    total_energy_faces::MaybeVector{Float64}
end

struct FluxVariables
    density_flux::Vector{Float64}
    momentum_flux::Vector{Float64}
    total_energy_flux::Vector{Float64}
end

include("../UI/initiate_UI.jl")
include("../UI/prompt_setup.jl")
include("../UI/SodShockUserInput.jl")

include("utility/plot_h5.jl")

include("grids/UniformAxis.jl")
include("grids/_1Dimension/Construct1DCartesian.jl")
include("grids/_1Dimension/Construct1DSpherical.jl")
include("grids/_2Dimension/Construct2DCartesian.jl")
include("grids/_2Dimension/Construct2DSpherical.jl")

include("grids/boundary_conditions/dispatch_boundary_condition.jl")
include("grids/boundary_conditions/outflow.jl")
include("grids/boundary_conditions/reflecting.jl")
include("grids/boundary_conditions/periodic.jl")

include("../examples/ShockTube/SodShockTube1D.jl")
include("../examples/ShockTube/RunShockTubeAnimate.jl")
include("../examples/ShockTube/RunShockTubeBenchmark.jl")
include("../examples/ShockTube/RunShockTubeDebug.jl")
include("../examples/ShockTube/RunShockTubeStandard.jl")
include("../examples/ShockTube/RunShockTubeVerbose.jl")

include("solvers/HYDRO/FDM/FTCS.jl")
include("solvers/HYDRO/FDM/LaxFriedrichs.jl")
include("solvers/HYDRO/FDM/Richtmyer.jl")


export run_simulation
export ConstructUniformAxis
export Construct1DCartesian
export Construct1DSpherical
export Construct2DCartesian
export Construct2DSpherical

export FTCS_Step
export LaxFriedrichs_Step
export Richtmyer_Step
#export guess_pressure
#export pressure_function
#export solve_star_region
#export sample_xi
export plot_snapshot
export animate_snapshots

function run_simulation()

 #   initiate_UI()
    UserInput = initiate_UI()
    
    if UserInput.primary_input.problem == :ShockTube
        if UserInput.primary_input.mode == :Animate
            run_Shock_Tube_Animate(UserInput)
        elseif UserInput.primary_input.mode == :Benchmark 
            run_Shock_Tube_Benchmark(UserInput)
        elseif UserInput.primary_input.mode == :Debug
            run_Shock_Tube_Debug(UserInput)
        elseif UserInput.primary_input.mode == :Standard
            run_Shock_Tube_Standard(UserInput) 
        elseif UserInput.primary_input.mode == :Verbose 
            run_Shock_Tube_Verbose(UserInput)
        end
        
    end
    
end

end
