function newton_divided_differences_logged(x, y; logfile::String = "Numerical_Methods/interpolation/debug/output/newton_debug.txt")
    open(logfile, "w") do io
        try
            write(io, "Newton Divided Differences Debug Log\n")
            write(io, "=====================================\n")
            write(io, "x = $(x)\n")
            write(io, "y = $(y)\n")

            n = length(x)
            F = zeros(n, n)
            F[:, 1] = y
            write(io, "Initial F = $(F)\n")

            for j in 2:n
                for i in 1:(n-j+1)
                    F[i, j] = (F[i+1, j-1] - F[i, j-1]) / (x[i+j-1] - x[i])
                end
                write(io, "F after column $j = $(F)\n")
            end

            write(io, "Final Coefficients = $(F[1, :])\n")
            return F[1, :], F

        catch e
            write(io, "\n--- ERROR ENCOUNTERED ---\n")
            write(io, "Error Type: $(typeof(e))\n")
            write(io, "Error Message: $(e.msg)\n")
            write(io, "Stacktrace:\n$(catch_backtrace())\n")
            rethrow()
        end
    end
end

function newton_interpolation_logged(x, y, x_eval; logfile::String = "Numerical_Methods/interpolation/debug/output/newton_debug.txt")
    open(logfile, "a") do io  # Append to same log
        try
            coeffs, _ = newton_divided_differences_logged(x, y; logfile=logfile)
            n = length(coeffs)
            write(io, "Evaluating at x_eval = $x_eval\n")
            write(io, "Coefficients = $(coeffs)\n")

            p = coeffs[n]
            for k in (n-1):-1:1
                p = coeffs[k] + (x_eval - x[k]) * p
            end
            write(io, "Interpolated value = $p\n")
            return p

        catch e
            write(io, "\n--- ERROR ENCOUNTERED DURING EVALUATION ---\n")
            write(io, "Error Type: $(typeof(e))\n")
            write(io, "Error Message: $(e.msg)\n")
            write(io, "Stacktrace:\n$(catch_backtrace())\n")
            rethrow()
        end
    end
end

function newton_interpolation_error_logged(x, y, x_eval; logfile::String = "Numerical_Methods/interpolation/debug/output/newton_debug.txt")
    open(logfile, "a") do io
        try
            _, F = newton_divided_differences_logged(x, y; logfile=logfile)
            n = length(x)
            max_derivative = abs(F[1, n])
            write(io, "Estimated max derivative = $max_derivative\n")

            error_term = 1.0
            for xi in x
                error_term *= (x_eval - xi)
            end

            error_bound = abs(error_term) * max_derivative / factorial(n)
            write(io, "Error term = $error_term\n")
            write(io, "Estimated error bound = $error_bound\n")

            return error_bound

        catch e
            write(io, "\n--- ERROR ENCOUNTERED DURING ERROR ESTIMATE ---\n")
            write(io, "Error Type: $(typeof(e))\n")
            write(io, "Error Message: $(e.msg)\n")
            write(io, "Stacktrace:\n$(catch_backtrace())\n")
            rethrow()
        end
    end
end
