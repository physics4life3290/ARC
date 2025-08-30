




function Codex_Trials(user_input::UserInput)

    init_bench_i, init_bench_f, prim_var_construct_i, prim_var_construct_f, evolution_bench_i, evolution_bench_f = nothing, nothing, nothing, nothing, nothing, nothing
    if :Benchmark ∈ user_input.Primary_Input.features
        init_bench_i, init_bench_f, prim_var_construct_i, prim_var_construct_f, evolution_bench_i, evolution_bench_f = Benchmark_Codex_Trials(user_input)
    end
    
    anim = Animation()
    h5_filename = user_input.Primary_Input.filename * ".h5"
    t_final = user_input.Solver_Input.t_final
    t = 0.0
    counter = 0
    _grid = nothing
    if user_input.Primary_Input.dimension == 1
        if user_input.Primary_Input.coordinate_system == :Cartesian
            _grid = Construct1DCartesian(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        elseif user_input.Primary_Input.coordinate_system == :Spherical || user_input.Primary_Input.coordinate_system == :Cylindrical
            _grid = Construct1DSpherical(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        end
    end
    
    println("The grid has been successfully created!")
    
    if user_input.Primary_Input.problem == :ShockTube

        if user_input.Primary_Input.dimension == 1
            # Problem Specific
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.states]
            densities = [state.density_centers for state in user_input.Secondary_Input.states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DShockTubePrimitives(_grid, user_input)

        end

    elseif user_input.Primary_Input.problem == :BlastWave

        if user_input.Primary_Input.dimension == 1
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.ambient_states]
            densities = [state.density_centers for state in user_input.Secondary_Input.ambient_states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DBlastWavePrimitives(_grid, user_input)
        end

    end

    U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, 
        W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (user_input.Secondary_Input.gamma - 1)) .+ W.velocity_centers.^2 ./ 2), 
        zeros(length(W.density_centers)), zeros(length(W.density_centers)), zeros(length(W.density_centers)))
    
    F = FluxVariables(zeros(length(U.density_centers)), 
            zeros(length(U.density_centers)),
            zeros(length(U.density_centers)))


    evolution_bench_i = time()
    while t < t_final
        counter += 1
        
        apply_boundary_conditions(user_input.Primary_Input.boundary_condition, U, _grid.coord1.zones, _grid.coord1.ghost_zones)
        
        if user_input.Primary_Input.solver == :FTCS
            FTCS_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :LaxFriedrichs
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :Richtmyer
            RichtmyerStep!(W, U, F, _grid, user_input,dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Secondary_Input.gamma, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :GodunovScheme
            GodunovStep!(W, U, F, user_input.Solver_Input.reconstruction, user_input.Solver_Input.limiter, user_input.Solver_Input.flattening, user_input.Solver_Input.steepening, user_input.Primary_Input.boundary_condition, user_input.Solver_Input.riemanntype, user_input.Secondary_Input.gamma, _grid.coord1.spacing, dt, user_input.Solver_Input.cfl, user_input.Primary_Input.mode, user_input.Primary_Input.features, _grid.coord1.total_zones, _grid.coord1.zones, _grid.coord1.ghost_zones)
        else 
            println("Defaulting to Lax until Scheme requested is supported...")
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        end

        source_terms = nothing
        if user_input.Primary_Input.coordinate_system == :cylindrical
            source_terms = (1/_grid.coord1.all_centers) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        elseif user_input.Primary_Input.coordinate_system == :spherical
            source_terms = (2/(_grid.coord1.all_centers)) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        end

        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2

        # This Clamping is for oscillatory schemes like FTCS and Richtmyer
        Threads.@threads for i in 1:_grid.coord1.total_zones
            @inbounds begin
                if W.pressure_centers[i] < 0.0
                    W.pressure_centers[i] = 0.0
                elseif W.density_centers[i] < 0.0
                    W.density_centers[i] = 0.0
                elseif W.internal_energy_centers[i] < 0.0
                    W.internal_energy_centers[i] = 0.0   
                end
            end
        end

        

        c = sqrt.(user_input.Secondary_Input.gamma .* W.pressure_centers ./ W.density_centers)
        dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / maximum(abs.(W.velocity_centers) .+ c)
        dt = min(dt, t_final - t)
        t += dt
   
        if counter % 10 == 0

            groupname = "step_$(counter)"
            println("Saving snapshot to $h5_filename in group $groupname")
            mode = isfile(h5_filename) ? "r+" : "w"
            h5open(h5_filename, mode) do file  # "a" = append mode
                grp = create_group(file, groupname)

                grp["x"] = _grid.coord1.all_centers

                grp["W/density"]      = W.density_centers
                grp["W/velocity"]     = W.velocity_centers
                grp["W/pressure"]     = W.pressure_centers
                grp["W/int_energy"]   = W.internal_energy_centers

                grp["U/density"]      = U.density_centers
                grp["U/momentum"]     = U.momentum_centers
                grp["U/total_energy"] = U.total_energy_centers

                grp["F/dens_flux"]    = F.density_flux
                grp["F/mome_flux"]    = F.momentum_flux
                grp["F/tot_ener_flux"]= F.total_energy_flux

                grp["time"] = t

            end
        
        end
        
    end

    function compute_norms(u_num::Vector, u_exact::Vector)
        diff = abs.(u_num .- u_exact)
        L1 = sum(diff)/length(diff)
        L2 = sqrt(sum(diff.^2)/length(diff))
        Linf = maximum(diff)
        return (L1=L1, L2=L2, Linf=Linf)
    end

    evolution_bench_f = time()
    if :Animate ∈ user_input.Primary_Input.features
        animate_snapshots(h5_filename, "W/density", savefile=user_input.Primary_Input.filename * ".gif")
    end
    println("The time it took to evolve the whole problem was: $(evolution_bench_f - evolution_bench_i) seconds...")
    dens_bench, vel_bench, press_bench = zeros(length(W.density_centers)), zeros(length(W.velocity_centers)), zeros(length(W.pressure_centers)) 

    if :Benchmark ∈ user_input.Primary_Input.features
        dens_bench .= W.density_centers
        vel_bench .= W.velocity_centers
        press_bench .= W.pressure_centers 

        if user_input.Primary_Input.problem == :ShockTube

            if user_input.Primary_Input.dimension == 1
                # Problem Specific
                pressures = [state.pressure_centers for state in user_input.Secondary_Input.states]
                densities = [state.density_centers for state in user_input.Secondary_Input.states]

                c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
                dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

                W = Construct1DShockTubePrimitives(_grid, user_input)

            end

        elseif user_input.Primary_Input.problem == :BlastWave

            if user_input.Primary_Input.dimension == 1
                pressures = [state.pressure_centers for state in user_input.Secondary_Input.ambient_states]
                densities = [state.density_centers for state in user_input.Secondary_Input.ambient_states]

                c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
                dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

                W = Construct1DBlastWavePrimitives(_grid, user_input)
            end

        end

        U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, 
            W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (user_input.Secondary_Input.gamma - 1)) .+ W.velocity_centers.^2 ./ 2), 
            zeros(length(W.density_centers)), zeros(length(W.density_centers)), zeros(length(W.density_centers)))
        
        F = FluxVariables(zeros(length(U.density_centers)), 
                zeros(length(U.density_centers)),
                zeros(length(U.density_centers)))

        function Riemann_Benchmark(x::Vector, xc::Vector, ρ::Vector, u::Vector, p::Vector, γ::Float64, t::Float64)
            N = length(x)
            M = length(xc)


            ρ_sol = similar(x)
            u_sol = similar(x)
            p_sol = similar(x)

            # Loop over spatial points where we want solution
            @inbounds for j in 1:M
                xj = x[j]

                # Find which cell this point is in
                cell = searchsortedlast(xc, xj)

                # Determine local Riemann problem (between cell-1 and cell)
                if cell == 0
                    # Left of domain: just take left state
                    ρ_sol[j], u_sol[j], p_sol[j] = ρ[1], u[1], p[1]
                    continue
                elseif cell == N
                    # Right of domain
                    ρ_sol[j], u_sol[j], p_sol[j] = ρ[N], u[N], p[N]
                    continue
                end

                # Interface position
                xif = 0.5 * (xc[cell] + xc[cell+1])

                # Left and right states
                ρL, uL, pL = ρ[cell], u[cell], p[cell]
                ρR, uR, pR = ρ[cell+1], u[cell+1], p[cell+1]

                # Solve Riemann problem for this interface
                pstar, ustar = solve_star_region(ρL, uL, pL, ρR, uR, pR, γ)

                ξ = (xj - xif) / t
                ρj, uj, pj = sample_exact(ξ, ρL, uL, pL, ρR, uR, pR, pstar, ustar, γ)

                ρ_sol[j] = ρj
                u_sol[j] = uj
                p_sol[j] = pj
            end

            return ρ_sol, u_sol, p_sol
        end

        while t < t_final

            counter += 1

            ρ, u, p = Riemann_Benchmark(_grid.coord1.centers, _grid.coord1.centers, W.density_centers, W.velocity_centers, W.pressure_centers, user_input.Secondary_Input.gamma, t)
            
            F.density_flux = ρ .* u
            F.momentum_flux = ρ .* u .^ 2 .+ p
            F.total_energy_flux = u .* (p ./ (user_input.Secondary_Input.gamma - 1) .+ 0.5 .* ρ .* u .^ 2 .+ p)

            source_terms = nothing
            if user_input.Primary_Input.coordinate_system == :cylindrical
                source_terms = (1/_grid.coord1.all_centers) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
            elseif user_input.Primary_Input.coordinate_system == :spherical
                source_terms = (2/(_grid.coord1.all_centers)) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
            end

            W.density_centers .= ρ
            W.velocity_centers .= u
            W.pressure_centers .= p
            W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2

            # This Clamping is for oscillatory schemes like FTCS and Richtmyer
            Threads.@threads for i in 1:_grid.coord1.total_zones
                @inbounds begin
                    if W.pressure_centers[i] < 0.0
                        W.pressure_centers[i] = 0.0
                    elseif W.density_centers[i] < 0.0
                        W.density_centers[i] = 0.0
                    elseif W.internal_energy_centers[i] < 0.0
                        W.internal_energy_centers[i] = 0.0   
                    end
                end
            end

            

            c = sqrt.(user_input.Secondary_Input.gamma .* W.pressure_centers ./ W.density_centers)
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / maximum(abs.(W.velocity_centers) .+ c)
            dt = min(dt, t_final - t)
            t += dt
        end
        dens_exact, vel_exact, press_exact = W.density_centers, W.velocity_centers, W.pressure_centers

        dens_norms = compute_norms(dens_bench, dens_exact)
        vel_norms = compute_norms(vel_bench, vel_exact)
        press_norms = compute_norms(press_bench, press_exact)

        println("Comparing the results from the $(user_input.Primary_Input.solver) solver to the exact Riemann solution: ")
        println("Density norms: ", dens_norms)
        println("Velocity norms: ", vel_norms)
        println("Pressure norms: ", press_norms)

    end

end

function Benchmark_Codex_Trials(user_input::UserInput)

    init_bench_i = time()
    t_final = user_input.Solver_Input.t_final
    t = 0.0
    _grid = nothing
    if user_input.Primary_Input.dimension == 1
        if user_input.Primary_Input.coordinate_system == :Cartesian
            _grid = Construct1DCartesian(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        elseif user_input.Primary_Input.coordinate_system == :Spherical || user_input.Primary_Input.coordinate_system == :Cylindrical
            _grid = Construct1DSpherical(user_input.Grid_Input.domain, user_input.Grid_Input.zones, user_input.Grid_Input.ghost_zones, user_input.Grid_Input.grid_center, "cm")
        end
    end
    prim_var_construct_i = time()
    if user_input.Primary_Input.problem == :ShockTube

        if user_input.Primary_Input.dimension == 1
            # Problem Specific
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.states]
            densities = [state.density_centers for state in user_input.Secondary_Input.states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DShockTubePrimitives(_grid, user_input)

        end

    elseif user_input.Primary_Input.problem == :BlastWave

        if user_input.Primary_Input.dimension == 1
            pressures = [state.pressure_centers for state in user_input.Secondary_Input.ambient_states]
            densities = [state.density_centers for state in user_input.Secondary_Input.ambient_states]

            c = sqrt(user_input.Secondary_Input.gamma * maximum(pressures) / maximum(densities))
            dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / c

            W = Construct1DBlastWavePrimitives(_grid, user_input)
        end

    end
    prim_var_construct_f = time()


    U = ConservativeVariables(W.density_centers, W.density_centers.* W.velocity_centers, 
        W.density_centers .* (W.internal_energy_centers .+ W.pressure_centers ./ (W.density_centers .* (user_input.Secondary_Input.gamma - 1)) .+ W.velocity_centers.^2 ./ 2), 
        zeros(length(W.density_centers)), zeros(length(W.density_centers)), zeros(length(W.density_centers)))
    
    F = FluxVariables(zeros(length(U.density_centers)), 
            zeros(length(U.density_centers)),
            zeros(length(U.density_centers)))

    init_bench_f = time()

    #if user_input.Primary_Input.mode == :Benchmark
    #    println("The time it took to create the initial primitive variables is: $(init_bench_f - init_bench_i) seconds...")
    #    println("The time it took to create the primitive variables is: $(prim_var_construct_f - prim_var_construct_i) seconds...")
    #end

    evolution_bench_i = time()
    while t < t_final
    
        apply_boundary_conditions(user_input.Primary_Input.boundary_condition, U, _grid.coord1.zones, _grid.coord1.ghost_zones)
        
        if user_input.Primary_Input.solver == :FTCS
            FTCS_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :LaxFriedrichs
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :Richtmyer
            RichtmyerStep!(W, U, F, _grid, user_input,dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Secondary_Input.gamma, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        elseif user_input.Primary_Input.solver == :GodunovScheme
            GodunovStep!(W, U, F, user_input.Solver_Input.reconstruction, user_input.Solver_Input.limiter, user_input.Solver_Input.flattening, user_input.Solver_Input.steepening, user_input.Primary_Input.boundary_condition, user_input.Solver_Input.riemanntype, user_input.Secondary_Input.gamma, _grid.coord1.spacing, dt, user_input.Solver_Input.cfl, user_input.Primary_Input.mode, user_input.Primary_Input.features, _grid.coord1.total_zones, _grid.coord1.zones, _grid.coord1.ghost_zones)
        else 
            println("Defaulting to Lax until Scheme requested is supported...")
            LaxFriedrichs_Step!(W, U, F, dt, _grid.coord1.ghost_zones, _grid.coord1.total_zones, _grid.coord1.spacing, user_input.Primary_Input.mode, user_input.Primary_Input.features, user_input.Solver_Input.cfl)
        end

        source_terms = nothing
        if user_input.Primary_Input.coordinate_system == :cylindrical
            source_terms = (1/_grid.coord1.all_centers) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        elseif user_input.Primary_Input.coordinate_system == :spherical
            source_terms = (2/(_grid.coord1.all_centers)) .* (F.density_flux, F.momentum_flux, F.total_energy_flux)
        end

        #if user_input.Primary_Input.mode == :Benchmark
        #    println("The time it took to apply the boundary conditions was: $(boundary_condition_bench_f - boundary_condition_bench_i) seconds...")
        #    println("The time it took to complete the solver step was: $(solver_step_bench_f - solver_step_bench_i) seconds...")
        #end

        W.density_centers .= U.density_centers
        W.velocity_centers .= U.momentum_centers ./ U.density_centers
        W.pressure_centers .= (user_input.Secondary_Input.gamma - 1) .* (U.total_energy_centers .- 0.5 .* U.density_centers .* W.velocity_centers .^ 2)
        W.internal_energy_centers .= U.total_energy_centers ./ U.density_centers .- 0.5 .* W.velocity_centers .^ 2

        # This Clamping is for oscillatory schemes like FTCS and Richtmyer
        Threads.@threads for i in 1:_grid.coord1.total_zones
            @inbounds begin
                if W.pressure_centers[i] < 0.0
                    W.pressure_centers[i] = 0.0
                elseif W.density_centers[i] < 0.0
                    W.density_centers[i] = 0.0
                elseif W.internal_energy_centers[i] < 0.0
                    W.internal_energy_centers[i] = 0.0   
                end
            end
        end

        c = sqrt.(user_input.Secondary_Input.gamma .* W.pressure_centers ./ W.density_centers)
        dt = user_input.Solver_Input.cfl * _grid.coord1.spacing / maximum(abs.(W.velocity_centers) .+ c)
        dt = min(dt, t_final - t)
        t += dt
        
    end
    evolution_bench_f = time()
    return init_bench_i, init_bench_f, prim_var_construct_i, prim_var_construct_f, evolution_bench_i, evolution_bench_f
end
