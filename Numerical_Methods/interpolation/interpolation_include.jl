

include("methods/cubic_spline.jl")
include("methods/divided_difference.jl")
include("methods/hermite_interp.jl")
include("methods/WENO_interp.jl")

include("debug/cubic_spline_debug.jl")
include("debug/divided_difference_debug.jl")
include("debug/hermite_interp_debug.jl")
include("debug/WENO_interp_debug.jl")

include("tests/interpolation_convergence_test.jl")
include("tests/cubic_spline_ptp_vs_all_points.jl")

include("utility/interp_input_check.jl")
include("utility/cubic_spline_dispatch.jl")
include("utility/div_diff_dispatch.jl")
include("utility/hermite_dispatch.jl")
include("utility/WENO_dispatch.jl")