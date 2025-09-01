




function get_regions(x, wall_positions, wall_iters)
    region_idx = 0
    if x <= wall_positions[1]
        region_idx = 1
    elseif x > wall_positions[end]
        region_idx = wall_iters + 1
    else
        for j in 1:wall_iters-1
            if wall_positions[j] < x <= wall_positions[j+1]
                region_idx = j + 1
                break
            end
        end
        if region_idx == 0
            error("Unable to assign region index to x = $x")
        end
    end
    return region_idx
end