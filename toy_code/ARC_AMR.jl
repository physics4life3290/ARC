using Printf
using Plots

# ----------------------
# Physics Constants
# ----------------------
const γ = 1.4

# ----------------------
# Block struct
# ----------------------
mutable struct Block
    # data is now a Vector of Vectors: [ [rho, rhou, E], ... ]
    data::Vector{Vector{Float64}}    
    dx::Float64
    level::Int
    children::Vector{Block}
    x0::Float64
    nx_internal::Int
end

internal_indices(b::Block) = 2:(b.nx_internal + 1)

# ----------------------
# Euler Physics Helpers
# ----------------------
function get_pressure(U)
    rho, rhou, E = U
    u = rhou / rho
    return (γ - 1.0) * (E - 0.5 * rho * u^2)
end

function get_primitive(U)
    rho = U[1]
    u = U[2] / rho
    p = get_pressure(U)
    return (rho, u, p)
end

function get_cons(rho, u, p)
    return [rho, rho * u, p / (γ - 1.0) + 0.5 * rho * u^2]
end

# HLL Riemann Solver for Fluxes
function hll_flux(UL, UR)
    rhoL, uL, pL = get_primitive(UL)
    rhoR, uR, pR = get_primitive(UR)
    
    aL = sqrt(γ * pL / rhoL)
    aR = sqrt(γ * pR / rhoR)
    
    # Simple wave speed estimates
    SL = min(uL - aL, uR - aR)
    SR = max(uL + aL, uR + aR)
    
    FL = [rhoL*uL, rhoL*uL^2 + pL, uL*(UL[3] + pL)]
    FR = [rhoR*uR, rhoR*uR^2 + pR, uR*(UR[3] + pR)]
    
    if SL >= 0
        return FL
    elseif SR <= 0
        return FR
    else
        return (SR .* FL .- SL .* FR .+ SL .* SR .* (UR .- UL)) ./ (SR - SL)
    end
end

# ----------------------
# Initialization & Refinement
# ----------------------
function create_block(nx, x0, dx, level, initial_U=[1.0, 0.0, 1.0])
    data = [copy(initial_U) for _ in 1:(nx + 2)]
    Block(data, dx, level, Block[], x0, nx)
end

function refine_block!(block::Block, threshold)
    empty!(block.children)
    nx = block.nx_internal
    idx = internal_indices(block)
    refine_mask = fill(false, nx)

    for i in 2:nx-1
        # 1. Variable access
        U = block.data[idx[i]]
        rho, p = U[1], get_pressure(U)
        
        # 2. Sensor A: Density Gradient (Relative)
        # Good for the contact discontinuity
        diff_rho = abs(block.data[idx[i+1]][1] - block.data[idx[i-1]][1]) / rho
        
        # 3. Sensor B: Pressure Gradient (Relative) 
        # Good for the expansion fan and shocks
        p_next, p_prev = get_pressure(block.data[idx[i+1]]), get_pressure(block.data[idx[i-1]])
        diff_p = abs(p_next - p_prev) / p

        # 4. Sensor C: Löhner Estimator (Density based)
        # Good for sharp shocks and avoiding over-refining smooth slopes
        num = abs(block.data[idx[i+1]][1] - 2*rho + block.data[idx[i-1]][1])
        den = abs(block.data[idx[i+1]][1] - rho) + abs(rho - block.data[idx[i-1]][1]) + 0.01*rho
        loehner = num / den

        # Combined "OR" logic with separate sensitivities
        if diff_rho > 0.1 || diff_p > 0.1 || loehner > threshold
            refine_mask[i] = true
            # Buffer zone: also refine neighbors to ensure the shock stays inside
            refine_mask[i-1] = true
            refine_mask[i+1] = true
        end
    end
    i = 1
    while i <= nx
        if refine_mask[i]
            start_i = i
            while i <= nx && refine_mask[i]; i += 1; end
            end_i = i - 1
            
            child_nx = (end_i - start_i + 1) * 2
            child_dx = block.dx / 2
            child_x0 = block.x0 + (start_i - 1) * block.dx
            child = create_block(child_nx, child_x0, child_dx, block.level + 1)
            
            for j in 1:child_nx
                parent_cell_idx = start_i + floor(Int, (j-1)/2)
                child.data[j+1] = copy(block.data[idx[parent_cell_idx]])
            end
            push!(block.children, child)
        else
            i += 1
        end
    end
end

# ----------------------
# Numerics
# ----------------------
function update_block!(block::Block, dt)
    nx = block.nx_internal
    idx = internal_indices(block)
    old_data = deepcopy(block.data)
    
    # Finite Volume Update
    for i in idx
        f_left = hll_flux(old_data[i-1], old_data[i])
        f_right = hll_flux(old_data[i], old_data[i+1])
        block.data[i] .-= (dt / block.dx) .* (f_right .- f_left)
    end
end

function restrict!(parent::Block)
    for child in parent.children
        p_idx = internal_indices(parent)
        c_idx = internal_indices(child)
        offset = round(Int, (child.x0 - parent.x0) / parent.dx)
        
        for j in 1:2:child.nx_internal
            parent_cell = offset + floor(Int, (j-1)/2) + 1
            parent.data[p_idx[parent_cell]] = 0.5 .* (child.data[c_idx[j]] .+ child.data[c_idx[j+1]])
        end
    end
end

function step_amr!(block::Block, dt)
    if isempty(block.children)
        update_block!(block, dt)
    else
        for _ in 1:2
            # Boundary: Child inherits from Parent (Simple)
            for child in block.children
                step_amr!(child, dt/2)
            end
        end
        restrict!(block)
    end
end

# ----------------------
# Execution
# ----------------------
function run_simulation()
    # 1. Setup Shock Tube Initial Conditions
    base = create_block(100, 0.0, 0.01, 0)
    for i in internal_indices(base)
        x = base.x0 + (i-1.5)*base.dx
        base.data[i] = (x < 0.5) ? get_cons(1.0, 0.0, 1.0) : get_cons(0.125, 0.0, 0.1)
    end

    dt = 0.001
    current_time = 0.0
    final_time = 0.5

    # 2. Capture the Animation object
    anim = @animate while current_time < final_time
        # Physics Step
        refine_block!(base, 0.02)
        step_amr!(base, dt)
        current_time += dt
        
        # Outflow Boundary Conditions
        base.data[1] = base.data[2]
        base.data[end] = base.data[end-1]

        # 3. Create the Frame with Annotation
        plt = plot(title="AMR Sod Shock Tube", ylims=(-0.1, 1.1), xlabel="x", ylabel="Density")
        plot_block_hierarchical!(plt, base)
        
        # Display current time on the frame
        annotate!(plt, 0.9, 1.0, text("t = $(round(current_time, digits=3))", :right, 10))
    end

    # 4. Save as MP4 instead of GIF
    # Note: Requires ffmpeg installed on your system
    mp4(anim, "shock_tube_simulation.mp4", fps = 20)
end



function plot_block_hierarchical!(plt, block::Block)
    idx = internal_indices(block)
    x = [block.x0 + (i-1.5)*block.dx for i in 1:block.nx_internal]
    y = [val[1] for val in block.data[idx]] # Plotting Density
    
    plot!(plt, x, y, lw=2, label="L$(block.level)", 
          marker=(block.level > 0 ? :circle : :none), ms=1.5)
    
    for child in block.children
        plot_block_hierarchical!(plt, child)
    end
end

run_simulation()
