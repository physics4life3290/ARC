




# Minmod slope limiter
function minmod(a, b)
    if a * b <= 0
        return 0.0
    else
        return sign(a) * min(abs(a), abs(b))
    end
end

@inline function minmod3(a, b, c)
    if (a > 0 && b > 0 && c > 0)
        return min(a, min(b, c))
    elseif (a < 0 && b < 0 && c < 0)
        return max(a, max(b, c))
    else
        return zero(a)
    end
end