module ARC

using CUDA
using HDF5
using Plots

include("utility/hydro_structs.jl")
include("grids/UniformAxis.jl")
include("grids/_1Dimension/Construct1DCartesian.jl")
include("grids/_1Dimension/Construct1DSpherical.jl")
include("grids/_2Dimension/Construct2DCartesian.jl")
include("grids/_2Dimension/Construct2DSpherical.jl")

include("solvers/HYDRO/FDM/FTCS.jl")
include("solvers/HYDRO/FDM/BTCS.jl")
include("solvers/HYDRO/FDM/LaxFriedrichs.jl")
include("solvers/HYDRO/FDM/Richtmyer.jl")

include("../examples/Linear1DAdvection.jl")
include("../examples/Linear2DAdvection.jl")
include("../examples/SodShockTube1D.jl")
include("../examples/SodShockTube2D.jl")
include("../examples/SedovBlastWave1D.jl")
include("../examples/SedovBlastWave2D.jl")
include("../examples/ShuOsher1D.jl")
include("../examples/ShuOsher2D.jl")
include("../examples/DoubleBlastWaveInteraction1D.jl")
include("../examples/DoubleBlastWaveInteraction2D.jl")
include("../examples/DoubleMachReflection2D.jl")


export run_simulation
export ConstructUniformAxis
export Construct1DCartesian
export Construct1DSpherical
export Construct2DCartesian
export Construct2DSpherical
export FTCS_Step
export BTCS_Step
export LaxFriedrichs_Step
export Richtmyer_Step

function run_simulation()

    _grid = Construct1DCartesian(2.0, 10, 3, 0.0, "cm")
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1
    γ = 4/3
    t_final = 0.2
    cfl = 0.3
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)
    dt = cfl * _grid.xcoord.spacing / max(cL, cR)

    W = Construct1DSodShockTubePrimitives(_grid.xcoord.total_zones, ρL, uL, pL, ρR, uR, pR, γ)
    U = Construct1DConservatives(W, γ)
    F = FluxVars(U.momentum, 
                 U.momentum .* W.velocity .+ W.pressure,
                 W.velocity .* (U.total_energy .+ W.pressure))

    RichtmyerStep(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing, γ)

    plot(_grid.xcoord.all_centers, W.density, label="Density", xlabel="Position (cm)", ylabel="Density")
    plot!(_grid.xcoord.all_centers, U.density, label="New Density")

end


end # module