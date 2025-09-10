





function enforce_monotonicity!(var::Vector{Float64}, varL::Vector{Float64}, varR::Vector{Float64})
    N = length(var)

    for i in 2:N-1
        v  = var[i]
        vL = varL[i]
        vR = varR[i]

        # -----------------------------
        # Case 1: local extremum
        # -----------------------------
        if (vR - v) * (v - vL) <= 0
            varL[i] = v
            varR[i] = v
            continue
        end

        # -----------------------------
        # Case 2: steep gradient correction
        # -----------------------------
        # left side adjustment
        if (vR - vL) * (vL - (3v - 2vR)) < 0
            varL[i] = 3v - 2vR
        end

        # right side adjustment
        if (vR - vL) * ((3v - 2vL) - vR) < 0
            varR[i] = 3v - 2vL
        end
    end
end
