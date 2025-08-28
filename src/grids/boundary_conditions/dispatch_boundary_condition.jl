


function apply_boundary_conditions(user_input, U, _grid)
    if user_input.Primary_Input.boundary_condition == :Reflecting
        apply_reflecting_boundaries!(U, _grid.coord1.ghost_zones, _grid.coord1.zones)
    elseif user_input.Primary_Input.boundary_condition == :Periodic
        apply_periodic_boundaries!(U, _grid.coord1.ghost_zones, _grid.coord1.zones)
    elseif user_input.Primary_Input.boundary_condition == :Outflow
        apply_outflow_boundaries!(U, _grid.coord1.ghost_zones, _grid.coord1.zones)
    else
        println("No boundary conditions...")
    end
end

function apply_boundary_conditions(boundary_condition::Symbol, U, zones::Int, ghost_zones::Int)
    if boundary_condition == :Reflecting
        apply_reflecting_boundaries!(U, ghost_zones, zones)
    elseif boundary_condition == :Periodic
        apply_periodic_boundaries!(U, ghost_zones, zones)
    elseif boundary_condition == :Outflow
        apply_outflow_boundaries!(U, ghost_zones, zones)
    else
        println("No boundary conditions...")
    end
end