


include("../methods/cubic_spline.jl")

# 1. Function the spline will ace (cubic polynomial)
f_ace(x) = x^3 #+ 2x^2 + x + 1

# 2. Function the spline will do okay on (smooth, non-polynomial)
f_okay(x) = sin(x)

# 3. Function that pushes the limits (discontinuous step function)
f_limit(x) = x < 0.5 ? 0.0 : 1.0


ind_var_exact = collect(range(0.0, 3.0, length=10000))
ind_var_interp = collect(range(0.1, 2.9, length = 10000))
dep_var_exact_ace = f_ace.(ind_var_exact)
dep_var_exact_okay = f_okay.(ind_var_exact)
dep_var_exact_limit = f_limit.(ind_var_exact)
ys = [dep_var_exact_ace, dep_var_exact_okay, dep_var_exact_limit]

interp_ys = []

for y in ys
    cs = cubic_spline(ind_var_exact, y)  # Create the cubic spline
    push!(interp_ys, cs.(ind_var_interp))  # Interpolate the values
end

using Plots

plot(ind_var_exact, dep_var_exact_ace ./ maximum(dep_var_exact_ace), label="Exact (Cubic)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
scatter!(ind_var_interp, interp_ys[1] ./ maximum(interp_ys[1]), label="Interpolated (Cubic)", linestyle=:dash)
plot!(ind_var_exact, dep_var_exact_okay ./ maximum(dep_var_exact_okay), label="Exact (Smooth)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
scatter!(ind_var_interp, interp_ys[2] ./ maximum(interp_ys[2]), label="Interpolated (Smooth)", linestyle=:dash)
plot!(ind_var_exact, dep_var_exact_limit ./ maximum(dep_var_exact_limit), label="Exact (Step)", title="Cubic Spline Interpolation", xlabel="x", ylabel="y")
scatter!(ind_var_interp, interp_ys[3] ./ maximum(interp_ys[3]), label="Interpolated (Step)", linestyle=:dash)
plot!(legend=:topright, size=(800, 600), grid=true)