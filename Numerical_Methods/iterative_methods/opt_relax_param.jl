

using LinearAlgebra

function opt_rel_param(D, L, U)
    
    T_gs = inv(D - L) .* U
    spectral_radius = maximum(eigen(T_gs).values)
    ω = 2 / (1 + sqrt(1-spectral_radius))
    
    return ω
end