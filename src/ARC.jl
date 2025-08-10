module ARC

using CUDA
using HDF5
using Plots
using ArgParse # For command line argument parsing access

const MaybeVector{T} = Union{Nothing, Vector{T}}

struct PrimitiveVariables 

    density_centers::Union{Vector{Float64}, Float64, Nothing}
    velocity_centers::Union{Vector{Float64}, Float64, Nothing}
    pressure_centers::Union{Vector{Float64}, Float64, Nothing}
    internal_energy_centers::MaybeVector{Float64}
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

include("../UI/UI_include.jl")

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

include("../examples/ShockTube/ShockTube1D.jl")
include("../examples/ShockTube/RunShockTubeAnimate.jl")
include("../examples/ShockTube/RunShockTubeBenchmark.jl")
include("../examples/ShockTube/RunShockTubeDebug.jl")
include("../examples/ShockTube/RunShockTubeStandard.jl")
include("../examples/ShockTube/RunShockTubeVerbose.jl")

include("solvers/HYDRO/FDM/FTCS.jl")
include("solvers/HYDRO/FDM/LaxFriedrichs.jl")
include("solvers/HYDRO/FDM/Richtmyer.jl")
include("solvers/HYDRO/FVM/GodunovStep.jl")
include("solvers/HYDRO/FVM/MUSCL.jl")
#include("solvers/HYDRO/FVM/PPM.jl")
include("solvers/HYDRO/ExactRiemannSolver.jl")
include("solvers/HYDRO/RiemannHLL.jl")
include("solvers/HYDRO/RiemannHLLC.jl")

export run_simulation
export ConstructUniformAxis
export Construct1DCartesian
export Construct1DSpherical
export Construct2DCartesian
export Construct2DSpherical

export FTCS_Step
export LaxFriedrichs_Step
export Richtmyer_Step
export ExactRiemannSolve!
export Godunov_Step!
export plot_snapshot
export animate_snapshots

function CalculateFlux!(W, U, F)
    F.density_flux .= U.momentum_centers 
    F.momentum_flux .= U.momentum_centers .* W.velocity_centers .+ W.pressure_centers
    F.total_energy_flux .= W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)
end



function run_simulation()

    #   initiate_UI()
    user_input = initiate_UI()
    println(user_input.Primary_Input)
    println(user_input.Grid_Input)
    println(user_input.Solver_Input)
    println(user_input.Secondary_Input)

    #=
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

    elseif UserInput.primary_input.problem == :BlastWave
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
    =#
    
end

end
