



function plot_snapshot(h5_filename::String, groupname::String; saveplot=false)
    h5open(h5_filename, "r") do file
        grp = file[groupname]

        x = grp["x"][]
        ρ = grp["W/density"][]
        v = grp["W/velocity"][]
        P = grp["W/pressure"][]
        e = grp["W/int_energy"][]
        t = grp["time"][]
        
        plt = plot(
            x, ρ ./ maximum(ρ), label="Density", lw=2, xlabel="x", ylabel="Normalized Value",
            title="Sod Shock Tube at t = $(round(t, digits=5)) s", legend=:topright
        )
        plot!(x, v ./ maximum(abs.(v)), label="Velocity", lw=2)
        plot!(x, P ./ maximum(P), label="Pressure", lw=2)
        plot!(x, e ./ maximum(e), label="Int. Energy", lw=2)

        if saveplot
            pngname = "plot_$(groupname).png"
            println("📸 Saving plot to $pngname")
            savefig(plt, pngname)
        else
            display(plt)
        end
    end
end

function animate_snapshots(h5_filename::String; gif_filename::String="animation.gif", fps::Int=30)
    groups = String[]

    # Collect all group names that start with "snapshot"
    h5open(h5_filename * ".h5", "r") do file
        for name in keys(file)
            if startswith(name, "snapshot")
                push!(groups, name)
            end
        end
    end

    sort!(groups)  # Ensure correct time order

    anim = @animate for groupname in groups
        h5open(h5_filename, "r") do file
            grp = file[groupname]

            x = grp["x"][]
            ρ = grp["W/density"][]
            v = grp["W/velocity"][]
            P = grp["W/pressure"][]
            e = grp["W/int_energy"][]
            t = grp["time"][]

            plot(
                x, ρ ./ maximum(ρ), label="Density", lw=2, xlabel="x", ylabel="Normalized Value",
                title="t = $(round(t, digits=5)) s", legend=:topright, ylim=(0, 1.1)
            )
            plot!(x, v ./ maximum(abs.(v)), label="Velocity", lw=2)
            plot!(x, P ./ maximum(P), label="Pressure", lw=2)
            plot!(x, e ./ maximum(e), label="Int. Energy", lw=2)
        end
    end

    println("Saving animation to $gif_filename")
    gif(anim, gif_filename, fps=fps)
end