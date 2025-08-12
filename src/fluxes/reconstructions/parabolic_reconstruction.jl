




function parabolic_reconstruct(var::AbstractVector; periodic::Bool=false)
    n = length(var)
    promT = promote_type(eltype(var), Float64)
    c1 = convert(promT, 7.0/12.0)
    c2 = convert(promT, -1.0/12.0)

    var_L = zeros(promT, n)
    var_R = zeros(promT, n)

    if periodic
        
        @inbounds for f in 2:n
            a = f-1                # left cell index
            b = f                  # right cell index
            im2 = (a - 2) % n + 1  # a-1 mod n
            ip1 = (b) % n + 1      # b+1 mod n
            face = c1*(promT(var[a]) + promT(var[b])) + c2*(promT(var[im2]) + promT(var[ip1]))
            var_L[f] = face
            var_R[f] = face
        end
        # fill face 1 (between cell n and 1) and face n+1 (same as face 1 for periodic)
        a = n; b = 1
        im2 = (a - 2) % n + 1
        ip1 = (b) % n + 1
        face01 = c1*(promT(var[a]) + promT(var[b])) + c2*(promT(var[im2]) + promT(var[ip1]))
        var_L[1] = face01
        var_R[1] = face01
        var_L[end] = face01
        var_R[end] = face01
    else
        
        for f in 2:n
            a = f-1
            b = f
            if a-1 < 1 || b+1 > n
                # stencil would go OOB -> fallback to simple average
                face = (promT(var[a]) + promT(var[b]))/2
            else
                face = c1*(promT(var[a]) + promT(var[b])) + c2*(promT(var[a-1]) + promT(var[b+1]))
            end
            var_L[f] = face
            var_R[f] = face
        end

        # boundary faces: set to cell value (no extrapolation)
        var_L[1] = promT(var[1]);   var_R[1] = promT(var[1])
        var_L[end] = promT(var[end]); var_R[end] = promT(var[end])
    end

    return var_L, var_R
end



function apply_ppm_limiter!(var::AbstractVector, var_L::AbstractVector, var_R::AbstractVector, user_input; periodic::Bool=false)
    n = length(var)
    #@assert length(var_L) == n+1 && length(var_R) == n+1 "face arrays must be length n+1"
    T = promote_type(eltype(var), eltype(var_L), Float64)

    # helper index functions for periodic wraps
    im1 = i -> periodic ? ((i - 2) % n) + 1 : (i - 1)
    ip1 = i -> periodic ? (i % n) + 1 : (i + 1)

    # 1) Clamp each face to [min(var[a],var[b]), max(...)], for interior faces
    for f in 2:n
        a = f-1; b = f
        lo = min(T(var[a]), T(var[b]))
        hi = max(T(var[a]), T(var[b]))
        var_L[f] = clamp(T(var_L[f]), lo, hi)
        var_R[f] = clamp(T(var_R[f]), lo, hi)
    end

    if periodic
        # clamp the wrap face (between cell n and 1)
        a, b = n, 1
        lo = min(T(var[a]), T(var[b]))
        hi = max(T(var[a]), T(var[b]))
        var_L[1] = clamp(T(var_L[1]), lo, hi)
        var_R[1] = clamp(T(var_R[1]), lo, hi)
        var_L[end] = var_L[1]
        var_R[end] = var_R[1]
    else
        # physical boundaries: no extrapolation
        var_L[1] = T(var[1])
        var_R[1] = T(var[1])
        var_L[end] = T(var[end])
        var_R[end] = T(var[end])
    end

    # 2) Per-cell check for envelope violation and slope limiting
    for i in 1:n-1
        # edge values coming from cell i
        q_left_from_i  = T(var_R[i])     # face i (left face) value from cell i
        q_right_from_i = T(var_L[i+1])   # face i+1 (right face) value from cell i

        # neighbors for envelope
        if periodic
            v_im1 = T(var[im1(i)])
            v_ip1 = T(var[ip1(i)])
        else
            v_im1 = (i == 1) ? T(var[i]) : T(var[i-1])
            v_ip1 = (i == n) ? T(var[i]) : T(var[i+1])
        end
        v_i = T(var[i])
        qmin = min(v_im1, min(v_i, v_ip1))
        qmax = max(v_im1, max(v_i, v_ip1))

        # Check if edges lie outside envelope
        if q_left_from_i < qmin || q_left_from_i > qmax || q_right_from_i < qmin || q_right_from_i > qmax
            # Compute differences for limiter
            if periodic
                Δplus  = T(var[ip1(i)]) - v_i
                Δminus = v_i - T(var[im1(i)])
                centered = (T(var[ip1(i)]) - T(var[im1(i)])) / 2
            else
                if i == 1
                    Δplus  = T(var[i+1]) - v_i
                    Δminus = zero(T)
                    centered = zero(T)
                elseif i == n
                    Δplus  = zero(T)
                    Δminus = v_i - T(var[i-1])
                    centered = zero(T)
                else
                    Δplus  = T(var[i+1]) - v_i
                    Δminus = v_i - T(var[i-1])
                    centered = (T(var[i+1]) - T(var[i-1])) / 2
                end
            end

            # Apply the chosen limiter
            s = zero(T)
            limiter = user_input.Solver_Input.limiter
            if limiter == :minmod
                s = minmod3(centered, 2 * Δplus, 2 * Δminus)
            elseif limiter == :superbee
                s = superbee(Δminus, Δplus)
            elseif limiter == :vanleer
                s = vanleer(Δminus, Δplus)
            else
                error("Unknown limiter: $limiter")
            end

            # Rebuild cell i's edges with slope-limited linear profile (MUSCL fallback)
            var_L[i+1] = v_i + 0.5 * s   # right edge from cell i
            var_R[i]   = v_i - 0.5 * s   # left edge from cell i
        end
    end

    return var_L, var_R
end
