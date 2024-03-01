using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings, Colors
include("../src/check_injectivity_radius.jl")

β = 0.5
Δρ = [1, 0.1, 0.01, 0.001, 0.0001]
iters = zeros(length(Δρ))
count = 1
coeff = solvef(π/2 + β, β)
for ρ ∈ (Δρ .+  min(π * coeff, min(sqrt(2β), 1) * π))
    global count, β
    bin, iters[count] =  check_radius(ρ, β, 4, 2, 10000000)
    count += 1
end
display(iters)
P = plot(framestyle=:box, legend=:topright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern", guidefontfamily="Computer Modern", titlefontfamily = "Computer Modern",
    legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13, titlefontsize = 18, margin = 0.3Plots.cm, xscale=:log, yscale=:log)

scatter!(Δρ, iters, markershape=:circle, markercolor =:black,  label=false, markersize=5, markerstrokewidth = 2)
plot!(Δρ, iters, linewidth = 2, label = false, color=:black, linestyle=:dash)
xlabel!(L"\rho-\hat{i}_\beta")
ylabel!("Iterations")
title!("Metric: " * L"\beta = "*string(β))
display(P)