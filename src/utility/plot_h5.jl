



function plot_snapshot(h5_filename::String, snapshot::Int, variable::String)
    # Construct the group name
    groupname = "step_$snapshot"
    
    # Open the HDF5 file
    h5open(h5_filename, "r") do file
        if haskey(file, groupname)
            grp = file[groupname]
            
            # Read x
            x = read(grp["x"])
            
            # Read the requested variable
            if !haskey(grp, variable)
                error("Variable '$variable' not found in group '$groupname'")
            end
            y = read(grp[variable])
            
            # Plot
            plot(x, y, xlabel="x", ylabel=variable, title="Snapshot $snapshot", lw=2)
        else
            error("Snapshot '$snapshot' not found in the file.")
        end
    end
end
using HDF5, Plots

function animate_snapshots(h5_filename::String, variables::Vector{String}; savefile::String="animation.mp4")
    h5open(h5_filename, "r") do file
        # Find all snapshot groups
        snapshots = sort([parse(Int, split(name, "_")[end]) 
                          for name in keys(file) if startswith(name, "step_")])

        # Create animation
        anim = @animate for snapshot in snapshots
            @inbounds begin
                println("Animating snapshot $snapshot out of $(length(snapshots))")
                grp = file["step_$snapshot"]

                # read x once per snapshot
                x = read(grp["x"])

                # Skip creating placeholder; only plot real variables
                first_plotted = false
                for var in variables
                    if !haskey(grp, var)
                        @warn "Variable '$var' not found in group 'step_$snapshot'"
                        continue
                    end
                    y = read(grp[var])
                    y_norm = y ./ maximum(abs.(y))    # normalize by own max

                    if !first_plotted
                        # First actual variable initializes the plot
                        plot(x, y_norm, label=var, xlabel="x", ylabel="Normalized value",
                             title="Snapshot $snapshot", lw=2)
                        first_plotted = true
                    else
                        # Subsequent variables are overplotted
                        plot!(x, y_norm, label=var)
                    end
                end
            end
        end

        # Save the animation as MP4
        mp4(anim, savefile, fps=30)
    end
end