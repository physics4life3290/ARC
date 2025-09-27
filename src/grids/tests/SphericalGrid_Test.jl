




include("../UniformAxis.jl")
include("../coordinate_systems/SphericalGrid.jl")

println("========== 1D Spherical Grid (Radial Shells) ==========")
grid1D = ConstructSphericalGrid(
    (10.0,),        # domain_lengths (r_max)
    (10,),          # zones
    3,              # ghost zones
    (0.0,)          # origin
)

println("Radial centers:      ", grid1D.axes[1].PhysicalCenters)
println("Radial spacing:      ", grid1D.axes[1].spacing)
println("Radial bounds:       ", grid1D.bounds[1])
println("Radial ghost bounds: ", grid1D.ghost_bounds[1])
println("Total zones:         ", grid1D.total_zones)
println("Origin:              ", grid1D.origin)
println("Dimension:           ", grid1D.dimension)

# ---------------------------------------------------------
println("\n========== 2D Spherical Grid (r–θ Plane) ==========")
grid2D = ConstructSphericalGrid(
    (10.0, Float64(π)),      # (r_max, θ_max)
    (10, 8),        # (radial zones, polar zones)
    3,
    (0.0, 0.0)
)

println("Radial centers:       ", grid2D.axes[1].PhysicalCenters)
println("Polar centers (θ):    ", grid2D.axes[2].PhysicalCenters)
println("Radial spacing:       ", grid2D.axes[1].spacing)
println("Polar spacing (θ):    ", grid2D.axes[2].spacing)
println("Bounds:               ", grid2D.bounds)
println("Ghost bounds:         ", grid2D.ghost_bounds)
println("Total zones:          ", grid2D.total_zones)
println("Origin:               ", grid2D.origin)
println("Dimension:            ", grid2D.dimension)

# ---------------------------------------------------------
println("\n========== 3D Spherical Grid (r–θ–φ) ==========")
grid3D = ConstructSphericalGrid(
    (10.0, Float64(π), Float64(2π)),  # (r_max, θ_max, φ_max)
    (10, 8, 16),    # (radial, polar, azimuthal)
    3,
    (0.0, 0.0, 0.0)
)

println("Radial centers:        ", grid3D.axes[1].PhysicalCenters)
println("Polar centers (θ):     ", grid3D.axes[2].PhysicalCenters)
println("Azimuthal centers (φ): ", grid3D.axes[3].PhysicalCenters)
println("Radial spacing:        ", grid3D.axes[1].spacing)
println("Polar spacing (θ):     ", grid3D.axes[2].spacing)
println("Azimuthal spacing (φ): ", grid3D.axes[3].spacing)
println("Bounds:                ", grid3D.bounds)
println("Ghost bounds:          ", grid3D.ghost_bounds)
println("Total zones:           ", grid3D.total_zones)
println("Origin:                ", grid3D.origin)
println("Dimension:             ", grid3D.dimension)
