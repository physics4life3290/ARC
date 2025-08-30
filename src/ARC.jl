module ARC

using CUDA
using HDF5
using Plots
using ArgParse # For command line argument parsing access
using Dates
using Printf

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
include("Codex_Trials.jl")

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
include("../examples/BlastWave1D.jl")

include("../logs/run_log.jl")
include("features/verbose/write_solver_input.jl")
include("features/verbose/write_solver_errors.jl")
include("features/verbose/write_solver_output.jl")

include("limiters/minmod.jl")
include("limiters/vanleer.jl")
include("limiters/superbee.jl")

include("fluxes/reconstructions/linear_reconstruction.jl")
include("fluxes/reconstructions/parabolic_reconstruction.jl")
include("fluxes/reconstructions/steepening.jl")
include("fluxes/reconstructions/flattening/flattening.jl")

include("solvers/HYDRO/FDM/FTCS.jl")
include("solvers/HYDRO/FDM/LaxFriedrichs.jl")
include("solvers/HYDRO/FDM/Richtmyer.jl")
include("solvers/HYDRO/FVM/GodunovScheme.jl")

include("fluxes/riemann_solvers/ExactRiemannSolver.jl")
include("fluxes/riemann_solvers/RiemannHLL.jl")
include("fluxes/riemann_solvers/RiemannHLLC.jl")

include("solvers/Interpolation/lagrangeinterp.jl")
include("solvers/Iterative/NewtonRaphson.jl")

export run_Codex_Trials
export ConstructUniformAxis
export Construct1DCartesian
export Construct1DSpherical
export Construct2DCartesian
export Construct2DSpherical

export FTCS_Step
export LaxFriedrichs_Step
export Richtmyer_Step
export Riemann_HLL
export Riemann_HLLC
export ExactRiemannSolve!
export Godunov_Step!
export MUSCL_Step!
export plot_snapshot
export animate_snapshots

function CalculateFlux!(W, U, F)
    F.density_flux .= U.momentum_centers 
    F.momentum_flux .= U.momentum_centers .* W.velocity_centers .+ W.pressure_centers
    F.total_energy_flux .= W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)
end

# This is where I need to begin and rehash the new modes, and their features.

function run_Codex_Trials()

    #   initiate_UI()
    user_input = initiate_UI()
    
    Codex_Trials(user_input)
    
end



end