




function get_blasts(P, x, γ,blasts)
    for blast in blasts
        dx = abs(x - blast.center)
        if dx <= blast.radius
            # Distribute blast energy as a simple top-hat pressure increase
            added_e = blast.energy / (2 * blast.radius)  # energy per unit length
            P += (γ - 1) * added_e
        end
    end
    return P
end