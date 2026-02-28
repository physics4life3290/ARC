# -----------------------------
# Cubic interpolation for 1D stencil
# -----------------------------
function cubic_interpolation_stencil(v::Vector{Float64}, dx::Vector{Float64})
    # Simple uniform-stencil cubic interpolation
    # v = [v_{i-1}, v_i, v_{i+1}, v_{i+2}]
    # dx = [dx_{i-1}, dx_i, dx_{i+1}, dx_{i+2}]
    # Returns: var_R, δvar_i, δvar_ip1

    # For uniform dx, standard PPM-style cubic reconstruction:
    var_R = (7*v[2] + 7*v[3] - v[1] - v[4])/12  # right interface of cell i
    δvar_i = (v[3] - v[1])/2                      # slope estimate for cell i
    δvar_ip1 = (v[4] - v[2])/2                    # slope estimate for cell i+1

    return var_R, δvar_i, δvar_ip1
end

function cubic_interpolation(var::Vector{Float64}, Δx::Vector{Float64}, i::Int)
    n = length(var)
    # build safe stencil indices
    im1 = max(i-1, 1)
    ip1 = min(i+1, n)
    ip2 = min(i+2, n)

    v_stencil = [var[im1], var[i], var[ip1], var[ip2]]
    dx_stencil = [
        Δx[im1],
        Δx[i],
        Δx[ip1],
        Δx[ip2]
    ]

    # Call your existing cubic interpolation stencil
    var_R, δvar_i, δvar_ip1 = cubic_interpolation_stencil(v_stencil, dx_stencil)
    return var_R, δvar_i, δvar_ip1
end

# -----------------------------
# Slope limiter (scalar version)
# -----------------------------
function slope_limiter(δvar, vi, vim1, vip1)
    if (vip1 - vi)*(vi - vim1) > 0
        return sign(δvar) * min(abs(δvar), 2abs(vi - vim1), 2abs(vip1 - vi))
    else
        return 0.0
    end
end

# -----------------------------
# Parabolic reconstruction in a cell
# -----------------------------
function parabolic_reconstruction(varL, varR, vi, ξ, ξL, ξR)
    α = (ξ - ξL) / (ξR - ξL)
    Δvar = varR - varL
    var6i = 6 * (vi - 0.5*(varL + varR))
    return varL + α*(Δvar + var6i*(1 - α))
end

# -----------------------------
# Reconstruct interfaces for N-cell array
# Returns left and right states
# -----------------------------
function reconstruct_interfaces(var::Vector{Float64}, coord::Vector{Float64})
    N = length(var)

    # Compute face coordinates
    ξf = zeros(N+1)
    Threads.@threads for i in 1:N-1
        @inbounds ξf[i+1] = 0.5*(coord[i] + coord[i+1])
    end
    # extrapolate outer faces
    ξf[1]   = coord[1]   - 0.5*(coord[2]-coord[1])
    ξf[end] = coord[end] + 0.5*(coord[end]-coord[end-1])

    # Compute cell widths
    Δξ = [ξf[i+1] - ξf[i] for i in 1:N]

    # Allocate interface arrays
    varL = zeros(N)
    varR = zeros(N)

    # Interior reconstruction
    for i in 2:N-2
        @inbounds begin
            # Cubic interpolation
            vR, δvar_i, δvar_ip1 = cubic_interpolation(var, Δξ, i)

            # Apply slope limiter
            δvar_i_lim = slope_limiter(δvar_i, var[i], var[i-1], var[i+1])

            # Assign right interface of cell i
            varR[i] = vR

            # Assign left interface of next cell (mirror)
            varL[i+1] = varR[i]
        end
    end

    # Boundaries: simple copy
    varL[1] = var[1]; varR[1] = var[1]
    varL[end] = var[end]; varR[end] = var[end]

    return varL, varR
end