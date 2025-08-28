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

export run_simulation
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

function run_simulation()

    #   initiate_UI()
    user_input = initiate_UI()

    init_bench_i = time()
    anim = Animation()
    Log = nothing
    h5_filename = user_input.Primary_Input.filename * ".h5"
        
    t_final = user_input.Solver_Input.t_final
    t = 0.0
    counter = 0
    _grid = nothing
    if user_input.Primary_Input.dimension == 1
        if user_input.Primary_Input.coordinate_system == :Cartesian
            _grid = Construct1DCartesian(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        elseif user_input.Primary_Input.coordinate_system == :Spherical || user_input.Primary_Input.coordinate_system == :Cylindrical
            _grid = Construct1DSpherical(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        end
    end
    println("The grid has been successfully created!")
    prim_var_construct_i = time()
    if user_input.Primary_Input.problem == :ShockTube

        if user_input.Primary_Input.dimension == 1
            # Problem Specific
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.states]
            densities = [state.density_centers for state in user_input.Secondary_Input.states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DShockTubePrimitives(_grid, user_input)

        end

    elseif user_input.Primary_Input.problem == :BlastWave

        if user_input.Primary_Input.dimension == 1
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.ambient_states]
            densities = [state.density_centers for state in user_input.Secondary_Input.ambient_states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DBlastWavePrimitives(_grid, user_input)
        end

    end
    prim_var_construct_f = time()


    U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, 
        W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (user_input.Secondary_Input.gamma - 1)) .+ W.velocity_centers.^2 ./ 2), 
        zeros(length(W.density_centers)), zeros(length(W.density_centers)), zeros(length(W.density_centers)))
    
    F = FluxVariables(zeros(length(U.density_centers)), 
            zeros(length(U.density_centers)),
            zeros(length(U.density_centers)))

    if user_input.Primary_Input.mode == :Verbose
        write_ShockTube_Log(user_input, Log) 
    end

    init_bench_f = time()

    if user_input.Primary_Input.mode == :Benchmark
        println("The time it took to create the initial primitive variables is: $(init_bench_f - init_bench_i) seconds...")
        println("The time it took to create the primitive variables is: $(prim_var_construct_f - prim_var_construct_i) seconds...")
    end

    evolution_bench_i = time()
    while t < t_final
        counter += 1
        boundary_condition_bench_i = time()
        apply_boundary_conditions(user_input.Primary_Input.boundary_condition, U, _grid.coord1.zones, _grid.coord1.ghost_zones)
        boundary_condition_bench_f = time()
        solver_step_bench_i = time()
        if user_input.Primary_Input.solver == :FTCS
            FTCS_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :LaxFriedrichs
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :Richtmyer
            RichtmyerStep!(W, U, F, _grid, user_input,dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Secondary_Input.gamma, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :GodunovScheme
            GodunovStep!(W, U, F, user_input.Solver_Input.reconstruction, user_input.Solver_Input.limiter, user_input.Solver_Input.flattening, user_input.Solver_Input.steepening, user_input.Primary_Input.boundary_condition, user_input.Solver_Input.riemanntype, user_input.Secondary_Input.gamma, _grid.coord1.spacing, dt, user_input.Solver_Input.cfl, user_input.Primary_Input.mode, user_input.Primary_Input.features, _grid.coord1.total_zones, _grid.coord1.zones, _grid.coord1.ghost_zones)
        else 
            println("Defaulting to Lax until Scheme requested is supported...")
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        end
        solver_step_bench_f = time()

        source_terms = nothing
        if user_input.Primary_Input.coordinate_system == :cylindrical
            source_terms = (1/_grid.coord1.all_centers) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        elseif user_input.Primary_Input.coordinate_system == :spherical
            source_terms = (2/(_grid.coord1.all_centers)) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        end

        if user_input.Primary_Input.mode == :Benchmark
            println("The time it took to apply the boundary conditions was: $(boundary_condition_bench_f - boundary_condition_bench_i) seconds...")
            println("The time it took to complete the solver step was: $(solver_step_bench_f - solver_step_bench_i) seconds...")
        end

        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2

        # This Clamping is for oscillatory schemes like FTCS and Richtmyer
        Threads.@threads for i in 1:_grid.coord1.total_zones
            @inbounds begin
                if W.pressure_centers[i] < 0.0
                    W.pressure_centers[i] = 0.0
                elseif W.density_centers[i] < 0.0
                    W.density_centers[i] = 0.0
                elseif W.internal_energy_centers[i] < 0.0
                    W.internal_energy_centers[i] = 0.0   
                end
            end
        end

        c = sqrt.(user_input.Secondary_Input.gamma .* W.pressure_centers ./ W.density_centers)
        dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / maximum(abs.(W.velocity_centers) .+ c)
        dt = min(dt, t_final - t)
        t += dt
   
        if counter % 10 == 0

            groupname = "step_$(counter)"
            println("Saving snapshot to $h5_filename in group $groupname")
            mode = isfile(h5_filename) ? "r+" : "w"
            h5open(h5_filename, mode) do file  # "a" = append mode
                grp = create_group(file, groupname)

                grp["x"] = _grid.coord1.all_centers

                grp["W/density"]      = W.density_centers
                grp["W/velocity"]     = W.velocity_centers
                grp["W/pressure"]     = W.pressure_centers
                grp["W/int_energy"]   = W.internal_energy_centers

                grp["U/density"]      = U.density_centers
                grp["U/momentum"]     = U.momentum_centers
                grp["U/total_energy"] = U.total_energy_centers

                grp["F/dens_flux"]    = F.density_flux
                grp["F/mome_flux"]    = F.momentum_flux
                grp["F/tot_ener_flux"]= F.total_energy_flux

                grp["time"] = t

            end
        
        end
        
    end
    evolution_bench_f = time()
    if :Animate ∈ user_input.Primary_Input.features
        animate_snapshots(h5_filename, "W/density", savefile=user_input.Primary_Input.filename * ".gif")
    end
    println("The time it took to evolve the whole problem was: $(evolution_bench_f - evolution_bench_i) seconds...")
end


end