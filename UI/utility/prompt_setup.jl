




# Helper functions for input
function prompt_choice_string(prompt::String, options::Vector{String}; allow_exit=true)
    while true
        println(prompt)
        for (i, opt) in pairs(options)
            println("$(i). $opt")
        end
        if allow_exit
            println("$(length(options)+1). Exit")
        end
        input = readline()
        choice = tryparse(Int, input)
        if choice !== nothing && 1 <= choice <= length(options) + (allow_exit ? 1 : 0)
            if allow_exit && choice == length(options)+1
                println("Exiting...")
                exit()
            end
            return choice
        elseif choice == "h" || choice == "help"
    
        else
            println("Invalid input. Please enter a number from the list.")
        end
    end
end

function prompt_number_string(prompt::String, T=Int)
    while true
        print(prompt * " ")
        input = readline()
        val = tryparse(T, input)
        if val !== nothing
            return val
        else
            println("Invalid input. Please enter a valid number.")
        end
    end
end

# Helper function to prompt user to choose from a list of Symbols
function prompt_choice(prompt::String, options::Vector{Symbol}; allow_exit=true)
    while true
        println(prompt)
        for (i, opt) in pairs(options)
            println("$(i). $opt")
        end
        if allow_exit
            println("$(length(options)+1). Exit")
        end

        print("> ")
        input = readline()
        choice = tryparse(Int, input)

        if choice !== nothing && 1 <= choice <= length(options) + (allow_exit ? 1 : 0)
            if allow_exit && choice == length(options) + 1
                println("Exiting...")
                exit()
            end
            return options[choice]
        else
            println("Invalid input. Please enter a number from the list.")
        end
    end
end

function prompt_choice_bool(prompt::String, options::Vector{Bool}; allow_exit=true)
    while true
        println(prompt)
        for (i, opt) in pairs(options)
            println("$(i). $opt")
        end
        if allow_exit
            println("$(length(options)+1). Exit")
        end

        print("> ")
        input = readline()
        choice = tryparse(Int, input)

        if choice !== nothing && 1 <= choice <= length(options) + (allow_exit ? 1 : 0)
            if allow_exit && choice == length(options) + 1
                println("Exiting...")
                exit()
            end
            return options[choice]
        else
            println("Invalid input. Please enter a number from the list.")
        end
    end
end

# Prompt for numeric input (unchanged)
function prompt_number(prompt::String, T=Int)
    while true
        print(prompt * " ")
        input = readline()
        val = tryparse(T, input)
        if val !== nothing
            return val
        else
            println("Invalid input. Please enter a valid number.")
        end
    end
end


function prompt_multiple_choices(prompt::String, options::Vector{Symbol}; allow_exit=true)
    while true
        println(prompt)
        for (i, opt) in pairs(options)
            println("$(i). $opt")
        end
        if allow_exit
            println("$(length(options)+1). Exit")
        end

        print("> ")
        input = readline()
        choices = try
            parse.(Int, split(input))
        catch
            println("Invalid input. Please enter space-separated numbers.")
            continue
        end

        # Handle Exit
        if allow_exit && (length(options)+1) in choices
            println("Exiting...")
            exit()
        end

        # Validate choices
        if all(c -> 1 <= c <= length(options), choices)
            unique_choices = unique(choices)
            return options[unique_choices]
        else
            println("Invalid input. Enter numbers between 1 and $(length(options)).")
        end
    end
end

function prompt_domain_length()
    while true
        println("Please input the domain length...")
        input = readline()
        domain_length = try
            parse(Float64, input)
        catch
            @error "Invalid input. Please enter a number."
            continue
        end

        if domain_length <= 0
            @error "Domain should be a positive number..."
            continue
        end

        return domain_length
    end
end

function prompt_num_states()
    while true
        println("Enter the number of states for the Shock Tube problem (at least 2):")
        input = readline()
        num_of_states = try
            parse(Int, input)
        catch
            @error "Invalid input. Please enter an integer."
            continue
        end

        if num_of_states < 2
            @error "There must be at least two states in a Shock Tube problem."
            continue
        end

        return num_of_states
    end
end

function prompt_adiabatic_constant()
    while true
        println("Please input an Adiabatic Constant (γ > 0, can be a fraction, e.g. 5/3):")
        input = strip(readline())
        value = try
            if occursin("/", input)
                nums = parse.(Float64, split(input, "/"))
                if length(nums) == 2 && nums[2] != 0.0
                    nums[1] / nums[2]
                else
                    @error "Invalid fraction. Please use 'numerator/denominator'."
                    continue
                end
            else
                parse(Float64, input)
            end
        catch
            @error "Invalid input. Please enter a number or fraction (e.g. 1.4 or 5/3)."
            continue
        end

        if value <= 0
            @error "The Adiabatic Constant must be greater than zero."
            continue
        end

        return value
    end
end

function prompt_num_blasts()
    while true
        println("How many blast regions do you want? (must be ≥ 1)")
        input = readline()
        num_blasts = try
            parse(Int, input)
        catch
            @error "Invalid input. Please enter an integer."
            continue
        end

        if num_blasts <= 0
            @error "Number of blast regions must be greater than zero."
            continue
        end

        return num_blasts
    end
end

function prompt_grid_min()
    while true
        println("Please input the minimum point of the grid (must be greater than 0):")
        input = strip(readline())
        value = try
            parse(Float64, input)
        catch
            @error "Invalid input. Please enter a number."
            continue
        end

        if value <= 0
            @error "The minimum point must be strictly greater than 0."
            continue
        end

        return value
    end
end

function prompt_zones()
    while true
        println("Please input the number of zones in the grid (must be > 0):")
        input = readline()
        zones = try
            parse(Int, input)
        catch
            @error "Invalid input. Please enter an integer."
            continue
        end

        if zones <= 0
            @error "Number of zones must be greater than zero."
            continue
        end

        return zones
    end
end

function prompt_cfl()
    while true
        println("Please input the Courant Condition (CFL > 0):")
        input = strip(readline())
        cfl = try
            parse(Float64, input)
        catch
            @error "Invalid input. Please enter a number."
            continue
        end

        if cfl <= 0
            @error "CFL must be greater than zero."
            continue
        end

        return cfl
    end
end

function prompt_t_final()
    while true
        println("Please input the final time (must be > 0):")
        input = strip(readline())
        t_final = try
            parse(Float64, input)
        catch
            @error "Invalid input. Please enter a number."
            continue
        end

        if t_final <= 0
            @error "Final time must be greater than zero."
            continue
        end

        return t_final
    end
end

