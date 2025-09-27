include("../UniformAxis.jl")
include("../coordinate_systems/CartesianGrid.jl")


grid1D = ConstructCartesianGrid((2.0,), (10,), 3, (0.0,))
println("=== 1D Grid ===")
println("Grid center: ", grid1D.origin[1])
println("Physical bounds: ", grid1D.bounds[1])
println("Ghost bounds: ", grid1D.ghost_bounds[1])
println("Total zones: ", grid1D.total_zones[1])



# 2D grid
grid2D = ConstructCartesianGrid(
    (10.0, 10.0),
    (10, 10),
    3,
    (0.0, 0.0)
)

println("\n=== 2D Grid ===")
println("Grid center: ", grid2D.origin)
println("Physical bounds: ", grid2D.bounds)
println("Ghost bounds: ", grid2D.ghost_bounds)
println("Total zones: ", grid2D.total_zones)


# 3D grid
grid3D = ConstructCartesianGrid(
    (2.0, 2.0, 2.0),
    (10, 10, 10),
    3,
    (0.0, 0.0, 0.0)
)

println("\n=== 3D Grid ===")
println("Grid center: ", grid3D.origin)
println("Physical bounds: ", grid3D.bounds)
println("Ghost bounds: ", grid3D.ghost_bounds)
println("Total zones: ", grid3D.total_zones)