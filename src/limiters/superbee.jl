



function superbee(a, b)
    if a * b <= 0
        return 0.0
    else
        s = sign(a)
        return s * max(0.0, min(2 * abs(a), abs(b)), min(abs(a), 2 * abs(b)))
    end
end