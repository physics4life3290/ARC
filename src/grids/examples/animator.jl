




using HDF5
using Plots

# Ensure the ffmpeg backend is available
gr()  # or pyplot(); ffmpeg will be used automatically by @animate

# Helper function to numerically sort step_* groups
function sort_groups_numerically(groupnames::Vector{String})
    steps = map(g -> parse(Int, split(g, "_")[2]), groupnames)
    sorted_indices = sortperm(steps)
    return groupnames[sorted_indices]
end

"""
animate_hdf5(h5_filename::String, dataset::String, variable::String, output_mp4::String; fps=15, lw=2)

Animates a variable stored in HDF5 snapshots saved under step_* groups.

# Arguments
- `h5_filename`: path to HDF5 file
- `dataset`: "W" or "U" (used to access variable like W/density)
- `variable`: the field to animate, e.g., "density", "velocity", "momentum"
- `output_mp4`: output filename for the MP4

# Optional keyword arguments
- `fps`: frames per second for the MP4 (default: 15)
- `lw`: line width for the plot (default: 2)
"""
function animate_hdf5(h5_filename::String, dataset::String, variable::String, output_mp4::String; fps=15, lw=2)
    # Open file in read-only mode
    h5file = h5open(h5_filename, "r")
    groupnames = sort_groups_numerically(collect(keys(h5file)))  # Proper numeric ordering

    # Get x-axis (assumes x stored in all snapshots identically)
    x = read(h5file[groupnames[1]]["x"])

    anim = @animate for groupname in groupnames
        grp = h5file[groupname]
        ykey = "$dataset/$variable"
        y = read(grp[ykey])
        time = read(grp["time"])

        plot(x, y,
            title = "$dataset/$variable at t = $(round(time, digits=3))",
            xlabel = "x",
            ylabel = variable,
            ylim = (0.0, 1.0),
            lw = lw,
            legend = false)
    end

    # Save the animation as MP4
    mp4(anim, output_mp4, fps = fps)

    close(h5file)
    println("Animation saved to $output_mp4")
end

animate_hdf5("test.h5", "W", "density", "density_evolution.mp4")
animate_hdf5("test.h5", "U", "momentum", "momentum_evolution.mp4")