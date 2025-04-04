





function hermite_interpolation(x_vals, f_vals, df_vals, x)
    n = length(x_vals)

    # Step 1: Expand the dataset
    X = repeat(x_vals, inner=2)  # Duplicate each x value
    F = repeat(f_vals, inner=2)  # Duplicate each f value

    # Step 2: Construct divided difference table
    H = zeros(2*n, 2*n)  # Hermite divided difference table
    H[:,1] = F           # First column is function values

    for i in 2:2:2n
        H[i,2] = df_vals[i÷2]  # First derivative entries
        if i > 2
            H[i-1,2] = (H[i-1,1] - H[i-2,1]) / (X[i-1] - X[i-2])
        end
    end

    for j in 3:2n
        for i in j:2n
            H[i,j] = (H[i,j-1] - H[i-1,j-1]) / (X[i] - X[i-j+1])
        end
    end

    # Step 3: Evaluate the Hermite polynomial using Newton’s form
    function newton_eval(x)
        result = H[1,1]
        product_term = 1.0

        for j in 2:2n
            product_term *= (x - X[j-1])
            result += H[j,j] * product_term
        end

        return result
    end

    return newton_eval(x)
end