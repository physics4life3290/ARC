



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

function animate_snapshots(h5_filename::String, variable::String; savefile::String="animation.gif")
    # Open HDF5 file
    h5open(h5_filename, "r") do file
        # Find all snapshot groups
        snapshots = sort([parse(Int, split(name, "_")[end]) 
                          for name in keys(file) if startswith(name, "step_")])

        # Create animation
        anim = @animate Threads.@threads for snapshot in snapshots
            @inbounds begin
                println("Animating snapshot $snapshot")
                grp = file["step_$snapshot"]
                
                if !haskey(grp, variable)
                    error("Variable '$variable' not found in group 'step_$snapshot'")
                end
                
                x = read(grp["x"])
                y = read(grp[variable])
                
                plot(x, y, xlabel="x", ylabel=variable, 
                    title="Snapshot $snapshot", lw=2)
            end
        end

        # Save the animation
        gif(anim, savefile, fps=60)
    end
end