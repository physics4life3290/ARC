



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