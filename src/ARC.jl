module ARC

using CUDA
using HDF5
using Plots
using ArgParse


include("../UI/initiate_UI.jl")
include("../UI/prompt_setup.jl")
include("../UI/choose_problem.jl")

include("utility/hydro_structs.jl")
include("utility/plot_h5.jl")

include("grids/UniformAxis.jl")
include("grids/_1Dimension/Construct1DCartesian.jl")
include("grids/_1Dimension/Construct1DSpherical.jl")
include("grids/_2Dimension/Construct2DCartesian.jl")
include("grids/_2Dimension/Construct2DSpherical.jl")

include("solvers/HYDRO/FDM/FTCS.jl")
include("solvers/HYDRO/FDM/BTCS.jl")
include("solvers/HYDRO/FDM/LaxFriedrichs.jl")
include("solvers/HYDRO/FDM/Richtmyer.jl")

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
export plot_snapshot

function ConservativetoPrimitive(U, W, γ)
    W.density = U.density
    W.velocity = U.momentum ./ U.density
    W.pressure = (γ - 1) * (U.total_energy - 0.5 * U.density .* W.velocity .^ 2)
    W.int_energy = U.total_energy ./ U.density - 0.5 * W.velocity .^ 2
end

function apply_reflecting_boundaries!(U::ConservativeVars, ng::Int, nx::Int)
    total = nx + 2ng

    # Left boundary
    for i in 1:ng
        U.density[ng - i + 1]        = U.density[ng + i]
        U.momentum[ng - i + 1]       = -U.momentum[ng + i]
        U.total_energy[ng - i + 1]   = U.total_energy[ng + i]
    end

    # Right boundary
    for i in 1:ng
        U.density[total - ng + i]       = U.density[total - ng - i + 1]
        U.momentum[total - ng + i]      = -U.momentum[total - ng - i + 1]
        U.total_energy[total - ng + i]  = U.total_energy[total - ng - i + 1]
    end
end

function apply_periodic_boundaries!(U::ConservativeVars, ng::Int, nx::Int)
    total = nx + 2ng

    # Fill left ghost zones (1:ng) with data from right physical zone (nx)
    for i in 1:ng
        U.density[ng - i + 1]        = U.density[nx + ng - i + 1]
        U.momentum[ng - i + 1]       = U.momentum[nx + ng - i + 1]
        U.total_energy[ng - i + 1]   = U.total_energy[nx + ng - i + 1]
    end

    # Fill right ghost zones (nx+ng+1 : total) with data from left physical zone (1:nx)
    for i in 1:ng
        U.density[nx + ng + i]       = U.density[ng + i]
        U.momentum[nx + ng + i]      = U.momentum[ng + i]
        U.total_energy[nx + ng + i]  = U.total_energy[ng + i]
    end
end

function apply_outflow_boundaries!(U::ConservativeVars, ng::Int, nx::Int)
    total = nx + 2ng

    # Left boundary: copy innermost left values into ghost zones
    for i in 1:ng
        U.density[ng - i + 1]        = U.density[ng + 1]
        U.momentum[ng - i + 1]       = U.momentum[ng + 1]
        U.total_energy[ng - i + 1]   = U.total_energy[ng + 1]
    end

    # Right boundary: copy innermost right values into ghost zones
    for i in 1:ng
        U.density[nx + ng + i]       = U.density[nx + ng]
        U.momentum[nx + ng + i]      = U.momentum[nx + ng]
        U.total_energy[nx + ng + i]  = U.total_energy[nx + ng]
    end
end

#=
function run_simulation()

    _grid = Construct1DCartesian(1.0, 500, 3, 0.0, "cm")
    ρL, uL, pL = 1.0, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1
    γ = 4/3
    t_final = 0.01
    t = 0.0
    cfl = 0.3
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)
    dt = cfl * _grid.xcoord.spacing / max(cL, cR)
    h5_filename = "SodShockTube1D.h5"
    
    W = Construct1DSodShockTubePrimitives(_grid.xcoord.total_zones, ρL, uL, pL, ρR, uR, pR, γ)
    #W = Construct1DShuOsherPrimitives(_grid.xcoord.all_centers, γ, 0.0)
    #W = Construct1DSedovBlastPrimitives(_grid.xcoord.all_centers, _grid.xcoord.spacing, 0.5, 1.0, 1.0, 1.0, γ)
    #W = Construct1DBlastWaveInteractionPrimitives(_grid.xcoord.all_centers, γ)
    
    U = Construct1DConservatives(W, γ)

    F = FluxVars(U.momentum, 
                U.momentum .* W.velocity .+ W.pressure,
                W.velocity .* (U.total_energy .+ W.pressure))

    counter = 0

    #anim = @animate while t < t_final 
    while t < t_final
        t += dt
        counter += 1
        #plot(title="Sod Shock Tube at t=$(t) s", xlabel="Position (cm)", ylabel="Density, Velocity, Pressure, Internal Energy", ylim=(0.0, 1.0), legend=:topright)

        println("Step ", counter, " Time: ", t)
        #plot!(_grid.xcoord.all_centers, W.density / maximum(W.density), label="Density")
        #plot!(_grid.xcoord.all_centers, W.velocity / maximum(abs.(W.velocity)), label="Velocity")
        #plot!(_grid.xcoord.all_centers, W.pressure / maximum(W.pressure), label="Pressure")
        #plot!(_grid.xcoord.all_centers, W.int_energy / maximum(W.int_energy), label="Internal Energy")

        #apply_reflecting_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
        #apply_periodic_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
        #apply_outflow_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
        #FTCS_Step!(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
        LaxFriedrichs_Step(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing)
        #RichtmyerStep!(U, F, dt, _grid.xcoord.ghost_zones, _grid.xcoord.total_zones, _grid.xcoord.spacing, γ)

        W.density .= U.density
        W.velocity .= U.momentum ./ U.density
        W.pressure .= (γ - 1) .* (U.total_energy .- 0.5 .* U.density .* W.velocity .^ 2)
        W.int_energy .= U.total_energy ./ U.density .- 0.5 .* W.velocity .^ 2
        

        # This Clamping is for oscillatory schemes like FTCS and Richtmyer
        for i in 1:_grid.xcoord.total_zones
            if W.pressure[i] < 0.0
                W.pressure[i] = 0.0
            elseif W.density[i] < 0.0
                W.density[i] = 0.0
            elseif W.int_energy[i] < 0.0
                W.int_energy[i] = 0.0   
            end
        end

        cs = sqrt.(γ .* W.pressure ./ W.density)
        wavespeed = abs.(W.velocity) .+ cs
        maxspeed = maximum(wavespeed)
        dt = cfl * _grid.xcoord.spacing / max(maxspeed, 1e-6)
        F.dens_flux .= U.momentum
        F.mome_flux .= U.momentum .* W.velocity .+ W.pressure
        F.tot_ener_flux .= W.velocity .* (U.total_energy .+ W.pressure)

        

        if counter % 1 == 0
            groupname = "step_$(counter)"
            println("Saving snapshot to $h5_filename in group $groupname")
            mode = isfile(h5_filename) ? "r+" : "w"
            h5open(h5_filename, mode) do file  # "a" = append mode
                grp = create_group(file, groupname)

                grp["x"] = _grid.xcoord.all_centers

                grp["W/density"]      = W.density
                grp["W/velocity"]     = W.velocity
                grp["W/pressure"]     = W.pressure
                grp["W/int_energy"]   = W.int_energy

                grp["U/density"]      = U.density
                grp["U/momentum"]     = U.momentum
                grp["U/total_energy"] = U.total_energy

                grp["F/dens_flux"]    = F.dens_flux
                grp["F/mome_flux"]    = F.mome_flux
                grp["F/tot_ener_flux"]= F.tot_ener_flux

                grp["time"] = t
            end
        end

    end

    #gif(anim, "Double_Blast_Wave_Interaction_Reflective.gif", fps=120)

end
=#

function run_simulation()

    initiate_UI()

end


end