module ARC

using CUDA
using HDF5
using Plots
using ArgParse
using Dates
using Printf


include("utility/hydro_structs.jl")
include("utility/exports.jl")
include("grids/grids_include.jl")
include("../UI/UI_include.jl")
include("Codex_Trials.jl")
include("utility/plot_h5.jl")
include("../examples/examples_include.jl")
include("../logs/run_log.jl")

include("solvers/solvers_include.jl")
include("fluxes/fluxes_include.jl")
include("features/features_include.jl")

function CalculateFlux!(W, U, F)
    F.density_flux .= U.momentum_centers 
    F.momentum_flux .= U.momentum_centers .* W.velocity_centers .+ W.pressure_centers
    F.total_energy_flux .= W.velocity_centers .* (U.total_energy_centers .+ W.pressure_centers)
end

function run_Codex_Trials()

    #   initiate_UI()
    user_input = initiate_UI()
    
    Codex_Trials(user_input)
    
end



end