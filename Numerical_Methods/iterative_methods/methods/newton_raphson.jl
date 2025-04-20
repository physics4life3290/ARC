
include("../interpolators/cubic_spline_interpolator.jl")

# I think I need to interpolate the function and its derivative.

function newton_raphson(f_interp, df_interp, x0; tol=1e-16, max_iter=100000, interp=true, legendre_order = 2)
    x = x0
    #iters = length(f)
    if interp == true
        for i in 1:max_iter
            #fx = f(x)
            #dfx = df(x)
            fx = f_interp(x)
            dfx = df_interp(x)
            #if abs(dfx) < tol
            #    error("Derivative too small, Newton-Raphson method may not converge.")
            #end
            x_new = x - fx / dfx
            if (abs(x_new - x)/abs(x)) < tol
                return x_new
            end
            x = x_new
        end
        error("Newton-Raphson method did not converge within the maximum number of iterations.")

    elseif interp != true

        for i in 1:max_iter
            f(x) = func(legendre_order, x)
            f_prime(x) = dfunc(legendre_order, x)
            
            # Check for convergence
            if abs(f) < tol
                return x
            end
            
            # Newton-Raphson update
            x_new = x - f / f_prime
            
            # Convergence check
            if abs(x_new - x) < tol
                return x_new
            end
            
            x = x_new
        end
        
        error("Newton-Raphson did not converge.")
    end

end


points = collect(1:1000)
func = [x^4 - 16 for x in 1:1000]
dfunc = [4*x^3 for x in 1:1000]
f_interp(var) = cubic_spline_interp(points, func, var)
df_interp(var) = cubic_spline_interp(points, dfunc, var)
root = newton_raphson(f_interp, df_interp, -1.0)
println("Root found: ", root)
