module ARC

using CUDA
using HDF5
using Plots

export run_simulation
export sarcastic_response

function run_simulation()
    println("ARC simulation running...")
    println("The sum of the first three primes is: $(2+3+5)")
    println("Did the Revise package work? I dont know, you tell me?")
end

function sarcastic_response()
    println("Naw, that shit dont work. That's some ole' bull shit!")
end

end # module