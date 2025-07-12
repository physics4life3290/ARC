module ARC

using CUDA
using HDF5
using Plots

include("grids/UniformAxis.jl")
include("grids/_1Dimension/Construct1DCartesian.jl")
include("grids/_1Dimension/Construct1DSpherical.jl")
include("grids/_2Dimension/Construct2DCartesian.jl")
include("grids/_2Dimension/Construct2DSpherical.jl")

export run_simulation
export ConstructUniformAxis
export Construct1DCartesian
export Construct1DSpherical
export Construct2DCartesian
export Construct2DSpherical

function run_simulation()
    
    cart1Dgrid = Construct1DCartesian(1.0, 100, 3, 0.5, "cm")

    t = 1.0
    CFL = 0.3
    # --- Initial condition ---
    function initial_condition(x)
        return exp.(-100 .* (x .- 0.5).^2)
    end

    function build1DLinearAdvection(cart1Dgrid, t, CFL)
        u = initial_condition(cart1Dgrid.xcoord.centers)
        u0 = copy(u)
        N = length(u)
        a = 1.0
        dt = CFL * cart1Dgrid.xcoord.spacing / abs(a)
        nt = Int(floor(t / dt))
        
        # --- Build BTCS matrix (periodic) ---
        r = a * dt / (2 * cart1Dgrid.xcoord.spacing)
        A = [i == j ? 1.0 : 0.0 for i in 1:N, j in 1:N]

        for i in 1:N
            A[i, mod1(i+1, N)] +=  r
            A[i, mod1(i-1, N)] += -r
        end
    end
    # --- Set up animation ---
    anim = @animate for n in 1:nt
        u = A \ u

        if n % 1 == 0 || n == 1 || n == nt
            plot(cart1Dgrid.xcoord.centers, u0, lw=2, label="Initial", xlabel="x", ylabel="u", legend=:topright)
            plot!(cart1Dgrid.xcoord.centers, u, lw=2, label="t = $(round(n*dt, digits=2))")
        end
    end

    # --- Save animation as GIF ---
    gif(anim, "btcs_advection.gif", fps=100)

end

end # module