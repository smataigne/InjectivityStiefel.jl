using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings, Colors
include("../src/check_injectivity_radius.jl")

#Choose n,p, n > p, and δ > 0 (the smaller δ, the longer the running time)
#Choose itermax
n = 5; p = 2
δ = 0.05
itermax = 1000

init = 0.35
dρ = [δ, 0]
βs = Array(init:0.05:1.5)

P = plot(framestyle=:box, legend=:bottomright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern", guidefontfamily="Computer Modern", titlefontfamily = "Computer Modern",
legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13, titlefontsize = 18, margin = 0.3Plots.cm)

for β ∈ βs
    global itermax
    coeff = solvef(π/2 + β, β)
    inj = (n - p > 1 ? min(π * coeff, min(sqrt(2β), 1) * π) : min(sqrt(2β), 1) * π)
    itermax = floor(Int, 6000 * exp(abs(β-0.5)/(β^2)))              #Customized itermax for quick run
    for ρ ∈ (dρ .+  inj)
        if check_radius(ρ, β, n, p, itermax)[1]
            if β == 0.1
                scatter!([β], [ρ], markershape=:circle, markercolor =:gray100,  label = "Stopped before iteration limit", markersize=4, markerstrokewidth = 2)
            else
                scatter!([β], [ρ], markershape=:circle, markercolor =:gray100, label = false, markersize=4, markerstrokewidth = 2)
            end
        else
            if β == 0.1
                scatter!([β], [ρ], markershape=:circle, color =:black, label = "Reached iteration limit", markersize=4)
            else
                scatter!([β], [ρ], markershape=:circle, color =:black, label = false, markersize=4)
            end
        end
    end
end

injs  = zeros(1000)
βs = LinRange(init, 1.5, 1000)
for i ∈ 1:1000
    β = βs[i]
    coeff = solvef(π/2 + β, β)
    injs[i] = (n - p > 1 ? min(π * coeff, min(sqrt(2β), 1) * π) : min(sqrt(2β), 1) * π)
end

plot!(βs, injs, label = L"\min\left\{\sqrt{2\beta}\pi, \pi, t_\beta^r \sqrt{2}\right\}", linewidth = 1, color =:black, linestyle=:solid)
xlabel!(L"\beta")
ylabel!(L"\rho")
display(P)
script_dir = @__DIR__
path = joinpath(script_dir, "../figures/plot_injectivity_beta_n_p_" *string(n) * "_" * string(p) * "_delta_"*string(δ) * ".pdf")
savefig(P, path)