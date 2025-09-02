




include("../limiters/minmod.jl")
include("../limiters/vanleer.jl")
include("../limiters/superbee.jl")

include("HYDRO/FDM/FTCS/FTCS.jl")
include("HYDRO/FDM/FTCS/Dispatch_FTCS.jl")
include("HYDRO/FDM/LaxFriedrichs/LaxFriedrichs.jl")
include("HYDRO/FDM/LaxFriedrichs/Dispatch_LaxFriedrichs.jl")
include("HYDRO/FDM/Richtmyer/Richtmyer.jl")
include("HYDRO/FDM/Richtmyer/Dispatch_Richtmyer.jl")
include("HYDRO/FVM/Godunov/GodunovScheme.jl")
include("HYDRO/FVM/Godunov/Dispatch_Godunov.jl")
include("Interpolation/lagrangeinterp.jl")
include("Iterative/NewtonRaphson.jl")
