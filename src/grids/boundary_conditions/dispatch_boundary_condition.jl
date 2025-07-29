




function apply_boundary_conditions(UserInput, U, _grid)
    if UserInput.primary_input.boundary_condition == :Reflecting
        apply_reflecting_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
    elseif UserInput.primary_input.boundary_condition == :Periodic
        apply_periodic_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
    elseif UserInput.primary_input.boundary_condition == :Outflow
        apply_outflow_boundaries!(U, _grid.xcoord.ghost_zones, _grid.xcoord.zones)
    end
end