using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings, Colors
include("check_injectivity_radius.jl")

δ = 0.07
dρ = [0, δ]
Δρ = []#[0.25,0.5, 0.75, 1]
βs = Array(0.1:0.05:1.5)

P = plot(framestyle=:box, legend=:bottomright,font="Computer Modern", tickfontfamily="Computer Modern",legendfont="Computer Modern", guidefontfamily="Computer Modern", titlefontfamily = "Computer Modern",
legendfontsize=12,yguidefontsize=17,xguidefontsize=17, xtickfontsize = 13, ytickfontsize=13, titlefontsize = 18, margin = 0.3Plots.cm)

for β ∈ βs
    coeff = solvef(π/2 + β, β)
    inj = min(π * coeff, min(sqrt(2β), 1) * π)
    for ρ ∈ (dρ .+  inj)
        if check_radius(ρ, β, 4, 2, 10000)[1]
            if β == 0.1
                scatter!([β], [ρ], markershape=:circle, markercolor =:gray100,  label="Stopped before iteration limit", markersize=5, markerstrokewidth = 2)
            else
                scatter!([β], [ρ], markershape=:circle, markercolor =:gray100, label=false, markersize=5, markerstrokewidth = 2)
            end
                #If check_radius succeeded for ρ, then it has succeeded for ρ + Δρ > ρ
            for ρ₂ ∈ (ρ .+ Δρ)
                scatter!([β], [ρ₂], markershape=:circle, color =:limegreen, label=false, markersize=5)
            end
        else
            if β == 0.1
                scatter!([β], [ρ], markershape=:circle, color =:black, label = "Reached iteration limit", markersize=5)
            else
                scatter!([β], [ρ], markershape=:circle, color =:black, label = false, markersize=5)
            end
            #If check_radius failed for ρ, then it has failed for ρ + Δρ < ρ
            for ρ₂ ∈ (ρ .- Δρ)
                scatter!([β], [ρ₂], markershape=skinnyx(0.2), color =:black, label = false, markersize=9)
            end
        end
    end
end

injs  = zeros(1000)
βs = LinRange(0.1, 1.5, 1000)
for i ∈ 1:1000
    β = βs[i]
    coeff = solvef(π/2 + β, β)
    injs[i] = min(π * coeff, min(sqrt(2β), 1) * π)
end

plot!(βs, injs, label = L"\min\left\{\sqrt{2\beta}\pi, \pi, t_\beta^r \sqrt{2}\right\}", linewidth = 1, color =:black, linestyle=:solid)
xlabel!(L"\beta")
ylabel!(L"\rho")
display(P)