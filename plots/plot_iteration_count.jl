using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings, Colors
include("../src/check_injectivity_radius.jl")

#Choose β>0, n>p>0
β = 0.8
n = 4 
p = 2
#Choose values of δ>0 to try.
N = 1
δ = exp2.(0:-1:-N)

iters = zeros(length(δ))
itermax  = 1
count = 1
coeff = solvef(π/2 + β, β)
for ρ ∈ (δ .+  min(π * coeff, min(sqrt(2β), 1) * π))
    global count, β
    bin, iters[count] =  check_radius(ρ, β, n, p, itermax)
    count += 1
end

P = plot(framestyle=:box, legend=:topright, font = "Computer Modern", tickfontfamily = "Computer Modern",legendfont = "Computer Modern", guidefontfamily="Computer Modern", titlefontfamily = "Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13, titlefontsize = 18, margin = 0.3Plots.cm, yscale=:log10, yticks = exp10.(0:Int(log10(itermax))))

scatter!(log2.(δ), iters, markershape=:circle, markercolor =:black,  label = false, markersize=5, markerstrokewidth = 2)
plot!(log2.(δ), iters, linewidth = 2, label = false, color=:black, linestyle=:dash)
xlabel!(L"\log_2(\rho-\hat{i}_\beta)")
ylabel!("Iterations")
title!("Metric: " * L"\beta = "*string(β))
display(P)

script_dir = @__DIR__
path = joinpath(script_dir, "../figures/Iterations_beta_" * string(β) * "_n_p_" *string(n) * "_" * string(p) * ".pdf")
#savefig(P, path)