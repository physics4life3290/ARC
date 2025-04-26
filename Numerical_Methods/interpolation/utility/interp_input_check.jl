




function prelim_interp_check(ind_var, dep_var)
    # Check if ind_var and dep_var are of the same type
    if typeof(ind_var) != typeof(dep_var)
        throw(ArgumentError("ind_var and dep_var must be of the same type"))
    end

    # Check if ind_var is sorted
    if !issorted(ind_var)
        throw(ArgumentError("ind_var must be sorted"))
    end

    # Check if ind_var and dep_var are not empty
    if isempty(ind_var) || isempty(dep_var)
        throw(ArgumentError("ind_var and dep_var must not be empty"))
    end
    
    # Check if the lengths of ind_var and dep_var are equal
    if length(ind_var) != length(dep_var)
        throw(ArgumentError("ind_var and dep_var must be the same length"))
    end
   
    # Check if ind_var and dep_var are of the same dimension
    if ndims(ind_var) != ndims(dep_var)
        throw(ArgumentError("ind_var and dep_var must be of the same dimension"))
    end

    # Check if ind_var and dep_var are of the same size
    if size(ind_var) != size(dep_var)
        throw(ArgumentError("ind_var and dep_var must be of the same size"))
    end
    
end