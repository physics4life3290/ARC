





function hermite_interpolation(x_vals, f_vals, df_vals, x)
    n = length(x_vals)
    @assert length(f_vals) == n && length(df_vals) == n

    # Expand x and f arrays
    X = repeat(x_vals, inner=2)
    H = zeros(2n, 2n)

    # Fill in function values and first divided differences
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

    # Fill in the rest of the divided difference table
    for j in 3:2n
        for i in j:2n
            H[i,j] = (H[i,j-1] - H[i-1,j-1]) / (X[i] - X[i-j+1])
        end
    end

    # Newton form evaluation
    function newton_eval(x_point)
        result = H[1,1]
        product_term = 1.0
        for j in 2:2n
            product_term *= (x_point - X[j-1])
            result += H[j,j] * product_term
        end
        return result
    end

    return newton_eval(x)
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