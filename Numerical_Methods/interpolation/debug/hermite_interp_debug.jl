




function hermite_interpolation_logged(x_vals, f_vals, df_vals, x_eval; logfile="Numerical_Methods/interpolation/debug/output/hermite_debug.txt")
    open(logfile, "w") do io
        try
            write(io, "Hermite Interpolation Debug Log\n")
            write(io, "===============================\n")
            write(io, "x_vals  = $x_vals\n")
            write(io, "f_vals  = $f_vals\n")
            write(io, "df_vals = $df_vals\n")
            write(io, "x_eval  = $x_eval\n\n")

            n = length(x_vals)
            @assert length(f_vals) == n && length(df_vals) == n

            X = repeat(x_vals, inner=2)
            H = zeros(2n, 2n)

            write(io, "Initial X = $X\n")
            write(io, "Initial H = $H\n\n")

            for i in 1:n
                H[2i-1, 1] = f_vals[i]
                H[2i,   1] = f_vals[i]
                H[2i,   2] = df_vals[i]
                if i == 1
                    H[2i-1, 2] = df_vals[i]
                else
                    H[2i-1, 2] = (H[2i-1,1] - H[2i-2,1]) / (X[2i-1] - X[2i-2])
                end
            end

            write(io, "H after filling first and second columns:\n$H\n\n")

            for j in 3:2n
                for i in j:2n
                    H[i,j] = (H[i,j-1] - H[i-1,j-1]) / (X[i] - X[i-j+1])
                end
            end

            write(io, "Final divided difference table H:\n$H\n\n")

            function newton_eval(x_point)
                result = H[1,1]
                product_term = 1.0
                for j in 2:2n
                    product_term *= (x_point - X[j-1])
                    result += H[j,j] * product_term
                end
                return result
            end

            interpolated = newton_eval(x_eval)
            write(io, "Interpolated value at x = $x_eval: $interpolated\n")
            return interpolated

        catch e
            write(io, "\n--- ERROR ENCOUNTERED ---\n")
            write(io, "Error Type: $(typeof(e))\n")
            write(io, "Error Message: $(e.msg)\n")
            write(io, "Stacktrace:\n$(catch_backtrace())\n")
            rethrow()
        end
    end
end


#=
x_vals = [1.0, 2.0, 3.0]
f(var) = var^5
df(var) = 5 * var^4
f_vals = f.(x_vals)
df_vals = df.(x_vals)
x = 1.5

println(hermite_interpolation(x_vals, f_vals, df_vals, x))
=#