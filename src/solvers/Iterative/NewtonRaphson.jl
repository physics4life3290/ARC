


#include("../Interpolation/lagrangeinterp.jl")

function newton_raphson(f, df, x0; tol=1e-10, maxiter=5000)
    x = x0
    #if typeof(f) == Function
    for k in 1:maxiter
        fx = f(x)
        dfx = df(x)
        if dfx == 0.0
            return x, false, k   # avoid divide by zero
        end
        dx = -fx / dfx
        x_new = x + dx
        # enforce positivity if needed
        if x_new <= 0
            x_new = 0.5 * x
        end
        x = x_new
        if abs(dx) < tol * (abs(x) + 1.0)
            return x, true, k
        end
    end
    return x, false, maxiter
end