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

function New_UI()
    modes = [:Standard, :Demonstrate, :Benchmark]
    problems = [:ShockTube, :BlastWave, :ShuOsher, :Custom]
    boundary_conditions = [:Reflective, :Periodic]
    grid_presets = [:Low, :Medium, :High, :Custom]
    adiabatic_const_presets = [:IdealGas, :Monatomic, :Diatomic, :Relativistic, :Custom]


    mode_prompt = ("Please select the mode you wish to run...")
    mode_choice = prompt_choice(mode_prompt, modes)

    problem_prompt = ("Please select the problem you wish to solve...")
    problem_choice = prompt_choice(problem_prompt, problems)

    println("Domain
    Default: [-1.0, 1.0]
    Enter new domain? (y/n)?")
    domain_choice = readline()
    
    if domain_choice == "y"
        println("Please enter the domain limits as two space-separated numbers (e.g. -0.5 0.5):")
        domain_input = readline()
        domain_limits = try
            parse.(Float64, split(domain_input))
        catch
            println("Invalid input. Using default domain [-1.0, 1.0].")
            [-1.0, 1.0]
        end
        if length(domain_limits) != 2 || domain_limits[1] >= domain_limits[2]
            println("Invalid domain limits. Using default domain [-1.0, 1.0].")
            domain_limits = [-1.0, 1.0]
        end
    else
        domain_limits = [-1.0, 1.0]
    end

    boundary_condition_prompt = ("Boundary Conditions")
    boundary_condition_choice = prompt_choice(boundary_condition_prompt, boundary_conditions)

    grid_preset_prompt = ("Grid Resolution Presets")


    if mode_choice != :Standard
        grid_preset_choice = prompt_choice(grid_preset_prompt, grid_presets[1:end-1])
    elseif mode_choice == :Standard
        grid_preset_choice = prompt_choice(grid_preset_prompt, grid_presets)
    end

    println("Final Time:
    Default: 1.0 seconds
    Enter new value? (y/n)")
    time_choice = readline()

    if time_choice == "y"
        println("Please enter the final time as a positive number (e.g. 0.5):")
        time_input = readline()
        final_time = try
            t = parse(Float64, time_input)
            if t <= 0
                println("Final time must be positive. Using default value of 1.0 seconds.")
                1.0
            else
                t
            end
        catch
            println("Invalid input. Using default final time of 1.0 seconds.")
            1.0
        end
    else
        final_time = 1.0
    end

    adiabatic_constant_prompt = ("Adiabatic Constant γ Presets:")
    if mode_choice != :Standard 
        adiabatic_constant_choice = prompt_choice(adiabatic_constant_prompt, adiabatic_const_presets[1:end-1])
    elseif mode_choice == :Standard
        adiabatic_constant_choice = prompt_choice(adiabatic_constant_prompt, adiabatic_const_presets)
    end
    
    println("Would you like to write a parameter file? (y/n)")
    file_choice = readline()

    println(mode_choice)
    println(problem_choice)
    println(domain_limits)
    println(boundary_condition_choice)
    println(grid_preset_choice)
    println(final_time)
    println(adiabatic_constant_choice)
end
export New_UI

end