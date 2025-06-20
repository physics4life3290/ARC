




function richtmyer_step!(U::AbstractMatrix, Δt, Δx)
    Unew = similar(U)
    half_step_U = similar(U)
    N = size(U,2)

    @views for i in 2:N-1
        half_step_U[:,i] .= 0.5*(U[:,i+1] + U[:,i-1]) .-
                            (Δt/(2Δx)) .* (flux(U[:,i+1]) .- flux(U[:,i-1]))
    end

    # Fix: assign edges for half_step_U to prevent undefined values
    @views half_step_U[:,1]   .= half_step_U[:,2]
    @views half_step_U[:,end] .= half_step_U[:,end-1]

    @views for i in 2:N-1
        Unew[:,i] .= U[:,i] .- (Δt/Δx) .* (flux(half_step_U[:,i+1]) .- flux(half_step_U[:,i-1]))
    end

    # simple outflow BC (zero-gradient)
    @views Unew[:,1]   .= Unew[:,2]
    @views Unew[:,end] .= Unew[:,end-1]

    return Unew
end