




include("../UniformAxis.jl")
include("../coordinate_systems/CylindricalGrid.jl")

include("../UniformAxis.jl")
include("../coordinate_systems/CylindricalGrid.jl")

# -------------------------------
# 1D cylindrical grid (radial only)
# -------------------------------
grid1D = ConstructCylindricalGrid(
    (10.0, ),  # (r_max, φ_max, z_max)
    (10, ),       # 10 radial zones, 1 azimuthal, 1 vertical
    3,                # ghost zones
    (0.0,)   # grid center
)

println("1D cylindrical grid (radial only):")
println("Radial centers: ", grid1D.axes[1].PhysicalCenters)
println("Radial spacing: ", grid1D.axes[1].spacing)
println("Radial bounds: ", grid1D.bounds[1])
println("Ghost bounds: ", grid1D.ghost_bounds[1])

# -------------------------------
# 2D cylindrical grid (r-φ plane)
# -------------------------------
grid2D = ConstructCylindricalGrid(
    (10.0, 2π),
    (10, 8),
    3,
    (0.0, 0.0)
)

println("\n2D cylindrical grid (r-φ plane):")
println("Radial centers: ", grid2D.axes[1].PhysicalCenters)
println("Azimuthal centers: ", grid2D.axes[2].PhysicalCenters)
println("Radial spacing: ", grid2D.axes[1].spacing)
println("Azimuthal spacing: ", grid2D.axes[2].spacing)
println("Bounds: ", grid2D.bounds)
println("Ghost bounds: ", grid2D.ghost_bounds)

# -------------------------------
# 3D cylindrical grid (r-φ-z)
# -------------------------------
grid3D = ConstructCylindricalGrid(
    (10.0, 2π, 5.0),
    (10, 8, 5),
    3,
    (0.0, 0.0, 0.0)
)

println("\n3D cylindrical grid (r-φ-z):")
println("Radial centers: ", grid3D.axes[1].PhysicalCenters)
println("Azimuthal centers: ", grid3D.axes[3].PhysicalCenters)
println("Axial centers: ", grid3D.axes[2].PhysicalCenters)
println("Spacings (r, φ, z): ", (grid3D.axes[1].spacing, grid3D.axes[3].spacing, grid3D.axes[2].spacing))
println("Bounds: ", grid3D.bounds)
println("Ghost bounds: ", grid3D.ghost_bounds)
