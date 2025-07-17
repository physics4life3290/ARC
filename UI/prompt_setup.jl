




# Helper functions for input
function prompt_choice(prompt::String, options::Vector{String}; allow_exit=true)
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
        else
            println("Invalid input. Please enter a number from the list.")
        end
    end
end

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