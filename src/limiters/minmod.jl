




# Minmod slope limiter
function minmod(a, b)
    if a * b <= 0
        return 0.0
    else
        return sign(a) * min(abs(a), abs(b))
    end
end

