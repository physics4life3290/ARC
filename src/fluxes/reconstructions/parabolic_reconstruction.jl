




function parabolic_reconstruct(var::AbstractArray)
    n = length(var)
    # Allocate left/right interface arrays
    var_L = zeros(n)
    var_R = zeros(n)

    # Parabolic reconstruction coefficients:
    # c0 = cell average (var[i])
    # c1 = slope term
    # c2 = curvature term

    Threads.@threads for i in 2:(n-1)
        @inbounds begin
            # Compute first differences
            d_im1 = var[i] - var[i-1]
            d_ip1 = var[i+1] - var[i]

            # Centered slope (first derivative)
            slope = 0.5 * (d_im1 + d_ip1)

            # Quadratic curvature (second derivative)
            curvature = 0.5 * (d_ip1 - d_im1)

            # Reconstruct left and right interface states
            # L state at i+1/2 comes from cell i
            var_R[i] = var[i] + 0.5 * slope + (1/6) * curvature

            # R state at i-1/2 comes from cell i
            var_L[i] = var[i] - 0.5 * slope + (1/6) * curvature
        end
    end

    # Boundaries: copy cell values
    var_L[1] = var[1]
    var_R[1] = var[1]
    var_L[end] = var[end]
    var_R[end] = var[end]

    return var_L, var_R
end



function apply_ppm_limiter!(var::AbstractVector, var_L::AbstractVector, var_R::AbstractVector, n::Int; user_input=nothing, periodic::Bool=false)

    @assert length(var_L) == n && length(var_R) == n "var_L/var_R must have length n (faces indexed 1..n)"
    @assert n >= 2 "PPM/MUSCL limiter needs at least 2 cells"

    T = promote_type(eltype(var), eltype(var_L), eltype(var_R), Float64)

    # cell index helpers
    im1 = i -> periodic ? ((i - 2) % n) + 1 : (i - 1)
    ip1 = i -> periodic ? (i % n) + 1       : (i + 1)

    # Face indexing convention used here:
    #   - Each cell i has a "left face index" fL(i) = i (interface between i-1 and i),
    #   - and a "right face index" fR(i) = periodic ? ip1(i) : (i < n ? i+1 : n).
    #   - In the periodic case, face f is between cells (f, ip1(f)).
    #   - In the non-periodic case, f=1 is the physical left boundary; f=n is the physical right boundary.

    # -----------------------------
    # 1) Clamp each face to local envelope [min(var[a],var[b]), max(...)]
    # -----------------------------
    if periodic
        # All faces are true interfaces
        Threads.@threads for f in 1:n
            @inbounds begin
                a = f
                b = ip1(f)
                lo = min(T(var[a]), T(var[b]))
                hi = max(T(var[a]), T(var[b]))
                var_L[f] = clamp(T(var_L[f]), lo, hi)
                var_R[f] = clamp(T(var_R[f]), lo, hi)
            end
        end
    else
        # Interior faces: 2..n-1 correspond to interfaces (i-1, i)
        Threads.@threads for f in 2:(n-1)
            @inbounds begin
                a = f - 1
                b = f
                lo = min(T(var[a]), T(var[b]))
                hi = max(T(var[a]), T(var[b]))
                var_L[f] = clamp(T(var_L[f]), lo, hi)
                var_R[f] = clamp(T(var_R[f]), lo, hi)
            end
        end
        # Physical boundaries: set to cell values (no extrapolation)
        var_L[1] = T(var[1]);   var_R[1] = T(var[1])   # left boundary
        var_L[n] = T(var[n]);   var_R[n] = T(var[n])   # right boundary
    end

    # -----------------------------
    # 2) Per-cell monotonicity check and slope-limited fallback (MUSCL)
    #    If either edge from cell i violates envelope, rebuild both edges from i.
    # -----------------------------

    if user_input !== nothing
        limiter = user_input.Solver_Input.limiter
    else
        limiter = lim 
    end

    # Helpers to safely fetch neighbor values in non-periodic case
    get_im1 = i -> periodic ? T(var[im1(i)]) : T(var[i == 1 ? 1 : i-1])
    get_ip1 = i -> periodic ? T(var[ip1(i)]) : T(var[i == n ? n : i+1])

    # Loop over cells; in non-periodic case we’ll avoid writing a non-existent right face
    istop = periodic ? n : n  # we still process i=n to possibly fix the left face
    Threads.@threads for i in 1:istop
        @inbounds begin
            v_i   = T(var[i])
            v_im1 = get_im1(i)
            v_ip1 = get_ip1(i)

            # Envelope across neighbors
            qmin = min(v_im1, min(v_i, v_ip1))
            qmax = max(v_im1, max(v_i, v_ip1))

            # Faces "owned" by cell i
            fL = i                                  # left face index
            fR = periodic ? ip1(i) : (i < n ? i+1 : n)  # right face index

            # Edge values coming from cell i
            q_left_from_i  = T(var_R[fL])     # left edge of cell i
            q_right_from_i = T(var_L[fR])     # right edge of cell i

            # If either edge lies outside the local envelope, rebuild using a slope-limited linear profile
            if (q_left_from_i < qmin || q_left_from_i > qmax ||
                q_right_from_i < qmin || q_right_from_i > qmax)

                # One-sided differences and centered difference
                if periodic
                    Δplus    = T(var[ip1(i)]) - v_i
                    Δminus   = v_i - T(var[im1(i)])
                    centered = (T(var[ip1(i)]) - T(var[im1(i)])) / 2
                else
                    if i == 1
                        Δplus    = T(var[i+1]) - v_i
                        Δminus   = zero(T)
                        centered = zero(T)
                    elseif i == n
                        Δplus    = zero(T)
                        Δminus   = v_i - T(var[i-1])
                        centered = zero(T)
                    else
                        Δplus    = T(var[i+1]) - v_i
                        Δminus   = v_i - T(var[i-1])
                        centered = (T(var[i+1]) - T(var[i-1])) / 2
                    end
                end

                # Choose limiter
                s::T = zero(T)
                if limiter == :minmod
                    s = minmod3(centered, 2*Δplus, 2*Δminus)
                elseif limiter == :superbee
                    s = superbee(Δminus, Δplus)
                elseif limiter == :vanleer
                    s = vanleer(Δminus, Δplus)
                else
                    error("Unknown limiter: $limiter")
                end

                # Rebuild edges for cell i (MUSCL fallback)
                # These are cell-i contributions to left and right faces
                var_R[fL] = v_i - 0.5*s
                if periodic || i < n
                    var_L[fR] = v_i + 0.5*s
                end

                # After rebuilding, clamp again to envelope to be safe
                var_R[fL] = clamp(var_R[fL], qmin, qmax)
                if periodic || i < n
                    var_L[fR] = clamp(var_L[fR], qmin, qmax)
                end
            end
        end
    end

    # -----------------------------
    # 3) Final boundary consistency
    # -----------------------------
    if periodic
        # Ensure wrap faces are clamped consistently between cells n and 1
        f = n # face between (n,1) in our periodic mapping
        lo = min(T(var[n]), T(var[1]))
        hi = max(T(var[n]), T(var[1]))
        var_L[f] = clamp(var_L[f], lo, hi)
        var_R[f] = clamp(var_R[f], lo, hi)
        # Face 1 already corresponds to (1,2) in periodic; nothing special to copy
    else
        # Hard-set physical boundaries (no extrapolation)
        var_L[1] = T(var[1]); var_R[1] = T(var[1])
        var_L[n] = T(var[n]); var_R[n] = T(var[n])
    end

    return var_L, var_R
end
