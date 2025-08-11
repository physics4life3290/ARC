




using Dates
using Printf


function write_ShockTube_Log(user_input, ShockTubeLog)
    write(ShockTubeLog, "
            #======================================#
            #   Sod Shock Tube / Riemann Problem   #
            #======================================#")

    write(ShockTubeLog, "
    
    User Defined Parameters:

                    Primary Input

        Problem               = $(user_input.Primary_Input.problem)
        Dimension             = $(user_input.Primary_Input.dimension)
        Coordinate System     = $(user_input.Primary_Input.coordinate_system)
        Solver                = $(user_input.Primary_Input.solver)
        Boundary Conditions   = $(user_input.Primary_Input.boundary_condition)
    
                    Secondary Input
                
        Domain                = $(user_input.Grid_Input.domain)
        Grid Center           = $(user_input.Grid_Input.grid_center)
        # of zones            = $(user_input.Grid_Input.zones)
        # of ghost zones      = $(user_input.Grid_Input.ghost_zones)
        Left Boundary         = $(user_input.Grid_Input.coord_min)
        Right Boundary        = $(user_input.Grid_Input.coord_max)
        # of States           = $(length(user_input.Secondary_Input.wall_positions) + 1)
        # of walls            = $(length(user_input.Secondary_Input.wall_positions))
        
                                ρ,   u,   P \n")

    if user_input.Primary_Input.problem == :ShockTube
        for i in 1:(length(user_input.Secondary_Input.wall_positions)+1)
            write(ShockTubeLog, "        State $i: $(user_input.Secondary_Input.states[i])\n")
        end
    elseif user_input.Primary_Input.problem == :BlastWave
        for i in 1:(length(user_input.Secondary_Input.wall_positions)+1)
            write(ShockTubeLog, "        State $i: $(user_input.Secondary_Input.parameters.states[i])\n")
        end
    end

    write(ShockTubeLog, "\n")
    
    for i in 1:(length(user_input.Secondary_Input.wall_positions))
        write(ShockTubeLog, "        Wall position $i: $(user_input.Secondary_Input.wall_positions[i])\n")
    end

    write(ShockTubeLog, "\n        Adiabatic Constant    = $(user_input.Secondary_Input.gamma)
        Courant #             = $(user_input.Solver_Input.cfl)
        Final Time            = $(user_input.Solver_Input.t_final)")
end

function write_State_ShockTube_log(_grid, W, U, F, ShockTubeLog)
    iters = length(_grid.xcoord.all_centers)
    write(ShockTubeLog, "\n                         Primitive Variables                                 Conservative Variables                         Flux Variables")
    write(ShockTubeLog, "\n   Radius        Density     Velocity     Pressure    Int. Energy        Density     Momentum    Tot.Energy         Density     Momentum    Tot.Energy\n")
    for i in 1:iters
        write(ShockTubeLog, @sprintf("%.6e %.6e %.6e %.6e %.6e     %.6e %.6e %.6e     %.6e %.6e %.6e\n",
        _grid.xcoord.all_centers[i],
        W.density_centers[i],
        W.velocity_centers[i],
        W.pressure_centers[i],
        W.internal_energy_centers[i], 
        U.density_centers[i],
        U.momentum_centers[i],
        U.total_energy_centers[i],
        F.density_flux[i],
        F.momentum_flux[i],
        F.total_energy_flux[i]))
    end
end