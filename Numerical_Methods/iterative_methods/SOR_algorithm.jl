



function sor(n::Int64, A::AbstractMatrix{Float64}, b::AbstractArray{Float64}, guess::AbstractArray{Float64}, ω::Float64, TOL::Float64, N::Int64)
    
    k = 1
    x_guess = copy(guess)  # Make a copy of the initial guess
    x_sol = zeros(Float64, n)
    temp = zeros(Float64, n)
    

    while k <= N
        #write(sor_log, "\n This is iteration $k...
        #The solution is...\n")
        @inbounds for i in 1:n
            temp[i] = dot(-A[i,1:i-1], x_sol[1:i-1]) + dot(-A[i,i+1:n], x_guess[i+1:n]) + b[i]
            dummy = x_sol[i] = ((1 - ω) * x_guess[i]) + (ω / A[i, i]) * temp[i]
            x_sol[1] = x_guess[1]
            x_sol[end] = x_guess[end]
            #write(sor_log, "$dummy \n")
        end

        _norm = norm(x_sol - x_guess)
        #write(sor_log, "The convergence check yields: $_norm...\n")

        if _norm < TOL
            #write(sor_log, "SOR Algorithm converged in $(k) iterations...\n")
            println("SOR Algorithm converged in $(k) iterations...\n")
            return x_sol
        end

        x_guess .= x_sol  # Update x_guess with the new solution
        k += 1
        
    end
    #write(main_sim_log, "The SOR Algorithm has encountered an error...\n")
    error("SOR did not converge in $N iterations...\n")
    #close(main_sim_log)
    #close(sor_log)
    #close(var_time_step_log)
end