




include("../make_grid.jl")
# ----------------------
# Example usage
# ----------------------
T = Float64
N = 2               # 1D
origin = (0.0, 0.0)     # note the comma!
domain = (1.0, 1.0)     # 1-element tuple
nx = (10, 10)          # must also be 1-element tuple
ghost_zones = 3

grid_u = make_grid(:uniform, origin, domain, nx, ghost_zones)

# Non-uniform grid
dx_array = (fill(0.1, nx[1]), fill(0.1, nx[2]))   # tuple of arrays
grid_nu = make_grid(:nonuniform, origin, domain, nx, ghost_zones, dx_array=dx_array)

# Adaptive grid
grid_ad = make_grid(:adaptive, origin, domain, nx, ghost_zones)

println(grid_u)
# ----------------------
# Parameters
# ----------------------
T = Float64
N = 1          # 1D
origin = (0.0,)
domain = (1.0,)
nx = (10,)     # 10 zones
ghost_zones = 2

# ----------------------
# Uniform grid
# ----------------------
grid_u = make_grid(:uniform, origin, domain, nx, ghost_zones)

println("Uniform Grid:")
println("Centers: ", grid_u.axes[1].centers)
println("Faces:   ", grid_u.axes[1].faces)
println("dx:      ", grid_u.axes[1].dx)
println("Bounds:  ", grid_u.axes[1].bounds)
println("Active zones: ", grid_u.active_zones[1])
println()

# ----------------------
# Non-uniform grid
# ----------------------
dx_array = (T[0.05,0.05,0.1,0.1,0.1,0.1,0.1,0.1,0.15,0.15],)
grid_nu = make_grid(:nonuniform, origin, domain, nx, ghost_zones, dx_array=dx_array)

println("Non-Uniform Grid:")
println("Centers: ", grid_nu.axes[1].centers)
println("Faces:   ", grid_nu.axes[1].faces)
println("dx:      ", grid_nu.axes[1].dx)
println("Bounds:  ", grid_nu.axes[1].bounds)
println("Active zones: ", grid_nu.active_zones[1])
println()

# ----------------------
# Adaptive grid (placeholder 1D block)
# ----------------------
grid_ad = make_grid(:adaptive, origin, domain, nx, ghost_zones, threshold=0.02)

println("Adaptive Grid:")
println("Root blocks: ", grid_ad.root_blocks)
println("Domain: ", grid_ad.domain)
println("Max level: ", grid_ad.max_level)
