




function apply_neumann_boundaries!(U::ConservativeVariables, ng::Int, nx::Int)
    total_centers = nx + 2ng

    # Left boundary (centers)
    left_val_idx = ng + 1
    @inbounds @simd for i in 1:ng
        dest = ng - i + 1
        U.density_centers[dest]      = U.density_centers[left_val_idx]
        U.momentum_centers[dest]     = U.momentum_centers[left_val_idx]
        U.total_energy_centers[dest] = U.total_energy_centers[left_val_idx]
    end

    # Right boundary (centers)
    right_val_idx = nx + ng
    @inbounds @simd for i in 1:ng
        dest = total_centers - ng + i
        U.density_centers[dest]      = U.density_centers[right_val_idx]
        U.momentum_centers[dest]     = U.momentum_centers[right_val_idx]
        U.total_energy_centers[dest] = U.total_energy_centers[right_val_idx]
    end

    # Faces (if they exist)
    if U.density_faces !== nothing
        total_faces = length(U.density_faces)
        
        left_face_idx = ng + 1
        right_face_idx = total_faces - ng
        
        # Left boundary (faces)
        @inbounds @simd for i in 1:ng
            dest = ng - i + 1
            U.density_faces[dest]      = U.density_faces[left_face_idx]
            U.momentum_faces[dest]     = U.momentum_faces[left_face_idx]
            U.total_energy_faces[dest] = U.total_energy_faces[left_face_idx]
        end

        # Right boundary (faces)
        @inbounds @simd for i in 1:ng
            dest = total_faces - ng + i
            U.density_faces[dest]      = U.density_faces[right_face_idx]
            U.momentum_faces[dest]     = U.momentum_faces[right_face_idx]
            U.total_energy_faces[dest] = U.total_energy_faces[right_face_idx]
        end
    end

    return nothing
end
