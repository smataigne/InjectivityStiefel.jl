using LinearAlgebra, SkewLinearAlgebra, Plots, LaTeXStrings
include("skewlog.jl")

@views function  dist_to_I_SOn(Q::AbstractMatrix)
    """
    Compute geodesic distance on SO(n)
    Input: Q, a matrix of SO(n)
    Output: Distance on SO(n)
    """
    return norm(skewlog(Q)) 
end

@views function givens(G::AbstractMatrix, ϕ::Real)
    """
    Builds Givens rotation of angle ϕ.
    """
    c = cos(ϕ); s = sin(ϕ)
    G[1, 1] = c; G[2, 1] = s
    G[1, 2] = -s; G[2, 2] = c
    return
end

"""
Equation to be solved for the bound on injectivity radius.
"""
f(t::Real, α::Real)  = sin(t) / t + (1 + 2α) * cos(t)
∇f(t::Real, α::Real) = (cos(t) * t - sin(t)) / (t * t) - (1 + 2α) * sin(t)

function solvef(t₀::Real, β::Real)
    """
    Newton's method to solve f(t, α) = 0
    """
    α = 1 / (2β) - 1
    δ = 1; tol = 1e-14
    itermax = 100; iter = 0
    tₙ  = t₀
    while abs(δ) > tol && iter < itermax
        δ = f(tₙ, α) / ∇f(tₙ, α)
        tₙ -= δ
        iter += 1
    end
    iter == itermax && return false
    return tₙ * sqrt(2) / π
end

@views function check_radius_canonical(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)
    """
    Returns true if ρ > inj. Stops at itermax otherwise and return false.
    """
    iszero(β - 1/2) || throw(ArgumentError("check_radius_canonical requires β = 1/2"))
    #Allocate memory
    Ω = zeros(n, n)
    A = zeros(p, p)
    B = zeros(n - p, p)
    temp = zeros(n, n - p)
    G = zeros(2, 2)
    iter = 0                              #Terminaison limit (set itermax = ∞ in theory)
    Nϕ = 1                                #Number of tries of rotation ϕ
    while iter < itermax
        #Instantiate unit random horizontal shooting direction
        A .= randn(p, p)                    
        skewhermitian!(A)
        B .= randn(n - p, p)
        nm = sqrt(β * norm(A)^2 + norm(B)^2)
        A ./= nm; B ./= nm
        Ω .= 0
        Ω[1:p, 1:p] .= A
        Ω[p+1:end, 1:p] .= B .* 2

        Q = exp(skewhermitian!(Ω) .* ρ)
        temp .= Q[:, p+1:end]

        for i ∈ 1:Nϕ                                  
            givens(G, rand() * 2π)
            mul!(Q[:, p+1:end], temp, G, 1, 0)     #Rotate last columns of Q
            if dist_to_I_SOn(Q) * sqrt(2) / 2 < ρ
                return true, iter + 1
            end
        end
        iter +=1
    end
    return false, iter
end

@views function check_radius(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)
    """
    Returns true if ρ > inj. Stops at itermax otherwise and return false.
    """
    β > 0 || throw(ArgumentError("check_radius_canonical requires β > 0"))
    iszero(β - 1/2) && return check_radius_canonical(ρ, β, n, p, itermax)
    #Allocate memory
    Ω = zeros(n, n)
    Ψ = zeros(p, p)
    A = zeros(p, p)
    B = zeros(n - p, p)
    tempQ1 = zeros(n, p)
    tempQ2 = zeros(n, n - p)
    tempV = zeros(p, p)
    G = zeros(2, 2)
    iter = 0                      #Terminaison limit (set itermax = ∞ in theory)
    Nϕ = 1                        #Number of tries of rotation ϕ
    while iter < itermax
        #Instantiate unit random horizontal shooting direction
        A .= randn(p, p)                    
        skewhermitian!(A)
        B .= randn(n - p, p)
        nm = sqrt(β * norm(A)^2 + norm(B)^2)
        A ./= nm; B ./= nm
        Ω .= 0
        Ω[1:p, 1:p] .=  A .* 2β
        Ω[p+1:end, 1:p] .= B .* 2
        skewhermitian!(Ω)
        Ψ .= A .* (2β - 1)
        Q = exp(skewhermitian!(Ω) .* ρ)
        V = exp(skewhermitian!(Ψ) .* ρ)

        tempQ1 .= Q[:, 1:p]
        tempQ2 .= Q[:, p+1:end]
        tempV  .= V

        for i ∈ 1:Nϕ    
            givens(G, rand() * 2π)  
            mul!(Q[:, 1:p], tempQ1, G, 1, 0)         #Rotate first columns of Q
            mul!(V, tempV, G, 1, 0)                  #Rotate V                                #
            givens(G, rand() * 2π)
            mul!(Q[:, p+1:end], tempQ2, G, 1, 0)     #Rotate last columns of Q 
            #Pseudo-Riemannian metric in total space if β > 0.5 
            Ω .=  skewlog(Q)
            Ψ .=  skewlog(V)                      
            d = sqrt(β * norm(Ω[1:p, 1:p] - Ψ)^2 + norm(Ω[p+1:end, 1:p])^2) 
            #d = sqrt((0.5 * dist_to_I_SOn(Q)^2 + β / (1 - 2β) * dist_to_I_SOn(V)^2))   
            if d < ρ
                return true , iter + 1
            end
        end
        iter +=1
    end
    return false , iter
end




