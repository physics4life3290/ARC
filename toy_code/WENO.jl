using Plots

# --- Constants ---
const γ = 1.4
const nx = 100
const x = range(-1.0, 1.0, length=nx)
const dx = step(x)
const t_final = 2.0  
const ε = 1e-10
const cfl = 0.3

get_p(ρ, m, E) = (γ - 1.0) * (E - 0.5 * m^2 / ρ)

function get_flux(U_vec)
    ρ, m, E = U_vec
    p = get_p(ρ, m, E)
    return [m, m^2/ρ + p, (E + p) * m / ρ]
end

# --- WENO5 with Reflective Ghost Cells ---
function weno5_reconstruct(v, var_type)
    n = length(v)
    
    # Create ghost cells (3 on each side for WENO5)
    # Left ghost cells: Mirror indices 3, 2, 1
    # Right ghost cells: Mirror indices n, n-1, n-2
    gL = [v[3], v[2], v[1]]
    gR = [v[n], v[n-1], v[n-2]]
    
    # If it's momentum (var_type == 2), we MUST negate the mirrored values
    if var_type == 2
        gL = -gL
        gR = -gR
    end
    
    # Combine: [Left Ghosts] + [Actual Data] + [Right Ghosts]
    ve = [gL; v; gR]
    
    v_out = zeros(n+1)
    # Offset is 3 because we added 3 ghost cells at the start
    for i in 4:n+4
        # Standard WENO5 reconstruction for interface i-1/2
        f1 = (1/3)*ve[i-3] - (7/6)*ve[i-2] + (11/6)*ve[i-1]
        f2 = -(1/6)*ve[i-2] + (5/6)*ve[i-1] + (1/3)*ve[i]
        f3 = (1/3)*ve[i-1] + (5/6)*ve[i] - (1/6)*ve[i+1]

        β1 = (13/12)*(ve[i-3]-2*ve[i-2]+ve[i-1])^2 + (1/4)*(ve[i-3]-4*ve[i-2]+3*ve[i-1])^2
        β2 = (13/12)*(ve[i-2]-2*ve[i-1]+ve[i])^2   + (1/4)*(ve[i-2]-ve[i])^2
        β3 = (13/12)*(ve[i-1]-2*ve[i]+ve[i+1])^2   + (1/4)*(3*ve[i-1]-4*ve[i]+ve[i+1])^2

        α1, α2, α3 = 0.1/(β1+ε)^2, 0.6/(β2+ε)^2, 0.3/(β3+ε)^2
        v_out[i-3] = (α1*f1 + α2*f2 + α3*f3) / (α1+α2+α3)
    end
    return v_out
end


function compute_rhs(U)
    rhs = [zeros(nx) for _ in 1:3]
    p = get_p.(U[1], U[2], U[3])
    u = U[2] ./ U[1]
    a = sqrt.(abs.(γ .* p ./ U[1]))
    α = maximum(abs.(u) .+ a)

    for j in 1:3
        f_cons = [get_flux([U[1][i], U[2][i], U[3][i]])[j] for i in 1:nx]
        fp, fm = 0.5 .* (f_cons .+ α .* U[j]), 0.5 .* (f_cons .- α .* U[j])
        
        # Use specific reflection logic for momentum (j=2)
        f_plus_inter  = weno5_reconstruct(fp, j)
        f_minus_inter = reverse(weno5_reconstruct(reverse(fm), j))
        
        flux_total = f_plus_inter .+ f_minus_inter
        for i in 1:nx
            # Flux at i+1/2 and i-1/2
            rhs[j][i] = -(flux_total[i+1] - flux_total[i]) / dx
        end
    end
    return rhs
end

# --- Simulation Initialization ---
ρ = [xi < 0.5 ? 1.0 : 0.125 for xi in x]
p = [xi < 0.5 ? 1.0 : 0.1 for xi in x]
U = [ρ, zeros(nx), p ./ (γ - 1.0)]

current_t, next_frame, frame_dt = 0.0, 0.0, 0.001

anim = @animate while current_t < t_final
    u, p = U[2]./U[1], get_p.(U[1], U[2], U[3])
    dt = cfl * dx / maximum(abs.(u) .+ sqrt.(abs.(γ .* p ./ U[1])))
    
    if current_t >= next_frame
        plot(x, U[1], title="Sod Shock: Reflective Walls (t=$(round(current_t, digits=3)))",
             ylims=(0, 1.1), lw=2, label="Density", xlabel="x")
        global next_frame += frame_dt
    end

    # SSP-RK3
    R1 = compute_rhs(U)
    U1 = [U[i] .+ dt .* R1[i] for i in 1:3]
    R2 = compute_rhs(U1)
    U2 = [0.75*U[i] + 0.25*U1[i] + 0.25*dt*R2[i] for i in 1:3]
    R3 = compute_rhs(U2)
    global U = [(1/3)*U[i] + (2/3)*U2[i] + (2/3)*dt*R3[i] for i in 1:3]
    global current_t += dt
end

mp4(anim, "sod_reflective_weno5.mp4", fps = 30)
