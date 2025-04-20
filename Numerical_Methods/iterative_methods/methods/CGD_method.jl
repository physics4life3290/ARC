using LinearAlgebra
using PyPlot

function conjugate_gradient(A, b, x0=nothing; tol=1e-5, max_iter=1000, boundaries=nothing)
    """
    Conjugate Gradient method for solving Ax = b with boundary conditions.

    Parameters:
    A : Matrix
        The matrix A in the system Ax = b
    b : Vector
        The right-hand side vector
    x0 : Vector, optional
        Initial guess for the solution (default is zero vector)
    tol : Float64, optional
        Tolerance for convergence (default is 1e-10)
    max_iter : Int, optional
        Maximum number of iterations (default is 1000)
    boundaries : Dict, optional
        Dictionary of boundary conditions {index: value} (default is nothing)

    Returns:
    x : Vector
        The solution vector
    """
    if x0 === nothing
        x0 = zeros(eltype(b), length(b))
    end

    if boundaries === nothing
        boundaries = Dict{Int, Float64}()
    end

    x = copy(x0)
    r = b - A * x
    for (index, value) in boundaries
        r[index] = 0.0  # Ensure residual respects boundary conditions
    end

    p = copy(r)
    rs_old = dot(r, r)

    for i in 1:max_iter
        Ap = A * p
        alpha = rs_old / dot(p, Ap)
        x = x .+ alpha .* p

        # Apply boundary conditions
        for (index, value) in boundaries
            x[index] = value
        end

        r = b - A * x
        for (index, value) in boundaries
            r[index] = 0.0  # Ensure residual respects boundary conditions
        end

        rs_new = dot(r, r)
        
        if sqrt(rs_new) < tol
            println("Converged in ", i, " iterations")
            break
        end

        p = r .+ (rs_new / rs_old) .* p
        for (index, value) in boundaries
            p[index] = 0.0  # Ensure direction respects boundary conditions
        end

        rs_old = rs_new
    end

    return x
end


