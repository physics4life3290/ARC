




function lax_friedrichs_step!(U::AbstractMatrix, Δt, Δx)
    Unew = similar(U)
    N = size(U,2)

    # interior update
    for i in 2:N-1
        Unew[:,i] .= 0.5*(U[:,i+1] + U[:,i-1]) .-
                      (Δt/(2Δx)).*(flux(U[:,i+1]) .- flux(U[:,i-1]))
    end

    # simple outflow BC (zero-gradient)
    Unew[:,1]   .= Unew[:,2]
    Unew[:,end] .= Unew[:,end-1]

    return Unew
end