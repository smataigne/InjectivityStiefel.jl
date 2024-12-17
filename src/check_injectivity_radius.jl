using LinearAlgebra, SkewLinearAlgebra, Distributions
include("skewlog.jl")

"""
dist_to_I_SOn(Q::AbstractMatrix)

Compute geodesic distance on SO(n)

Input: Q, a matrix of SO(n)
Output: Distance on SO(n)
"""
@views function  dist_to_I_SOn(Q::AbstractMatrix)
    return norm(skewlog(Q)) 
end

"""
givens(G::AbstractMatrix, ϕ::Real)

Builds Givens rotation of angle ϕ.
"""
@views function givens(G::AbstractMatrix, ϕ::Real)
   
    c = cos(ϕ); s = sin(ϕ)
    G[1, 1] = c; G[2, 1] = s
    G[1, 2] = -s; G[2, 2] = c
    return
end

"""
f(t::Real, α::Real)

f(t) = 0 is the equation to be solved for the bound on injectivity radius with the α-metric.
"""
f(t::Real, α::Real)  = sin(t) / t + (1 + 2α) * cos(t)
∇f(t::Real, α::Real) = (cos(t) * t - sin(t)) / (t * t) - (1 + 2α) * sin(t)

"""
solvef(t₀::Real, β::Real)

Newton's method to solve f(t, α) = 0
"""
function solvef(t₀::Real, β::Real)
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

"""
check_radius_canonical(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)

Returns true if ρ ter than the injectivity radius on the Stiefel manifold St(n,p) endowed with the β-metric.

Stops at itermax otherwise and return false.
"""
@views function check_radius_canonical(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)
    iszero(β - 1/2) || throw(ArgumentError("check_radius_canonical requires β = 1/2"))
    #If canonical metric and n - p = 1, the search on the fiber makes no sense.
    iszero(n - p - 1) && return ρ > π, 1
    #Allocate memory
    Ω = zeros(n, n)
    A = zeros(p, p)
    B = zeros(n - p, p)
    temp = zeros(n, n - p)
    G = zeros(n - p, n - p)
    iter = 0                              #Terminaison limit (set itermax = ∞ in theory)
    Nϕ = 5                                #Number of initial guesses on the fiber
    α = 0.5                               #Stepsize for gradient method
    ∇tol = 1e-3; ∇max = 40                #Gradient descent parameters
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
            if n - p == 2              
                givens(G, randn() * 2π)
            else
                G = exp(skewhermitian!(randn(n - p, n - p)))
            end
            mul!(Q[:, p+1:end], temp, G, 1, 0)     #Rotate last columns of Q
            #Gradient descent over squared distance on the fiber
            ∇iter = 0
            ∇Q = skewlog(Q)[p+1:end, p+1:end]
            while norm(∇Q) > ∇tol && ∇iter < ∇max
                temp .= Q[:, p+1:end]
                ∇Q = skewlog(Q)[p+1:end, p+1:end]
                G = exp(-α * ∇Q)
                mul!(Q[:, p+1:end], temp, G, 1, 0)     #Rotate last columns of Q
                ∇iter += 1
            end
            if dist_to_I_SOn(Q) * sqrt(2) / 2 < ρ
                print("ρ is greater than the injectivity radius \n")
                return true, iter + 1
            end
        end
        iter +=1
    end
    print("Trial budget exceeded. Either ρ is smaller than the injectivity radius; 
    or it is larger and the algorithm has not been able to find a closer cut time, 
    a likely event when ρ is only slightly larger than the injectivity radius. \n")
    return false, iter
end

@views function choose_logarithm!(V::AbstractVector, K::AbstractVector)
    for i ∈ 1:2:length(V)
        V[i] .+= K[(i+1)÷2] * 1im * 2π
        V[i + 1] .-= K[(i+1)÷2] * 1im * 2π
    end
end

"""
check_radius(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)

Returns true if ρ ter than the injectivity radius on the Stiefel manifold St(n,p) endowed with the β-metric.
     
Stops at itermax otherwise and return false.
"""
@views function check_radius(ρ::Real, β::Real, n::Integer, p::Integer, itermax::Integer)
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
    G1 = zeros(p, p)
    G2 = zeros(n - p, n - p)
    iter = 0                       #Terminaison limit (set itermax = ∞ in theory)
    Nϕ = 1                         #Number of initial guesses on the fiber
    α = 0.5                        #Stepsize for gradient method
    ∇tol = 1e-1; ∇max = 100        #Gradient descent parameters
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
            if p == 2    
                givens(G1, rand() * 2π)  
            else
                G1 = exp(skewhermitian!(randn(p, p) * 2π))
            end
            mul!(Q[:, 1:p], tempQ1, G1, 1, 0)         #Rotate first columns of Q
            mul!(V, tempV, G1, 1, 0)                  #Rotate V    
            if n - p == 2              
                givens(G2, rand() * 2π)
            else
                G2 = exp(skewhermitian!(randn(n - p, n - p)))
            end
            mul!(Q[:, p+1:end], tempQ2, G2, 1, 0)     #Rotate last columns of Q 
            #Pseudo-Riemannian metric in total space if β > 0.5 
            #Gradient descent over squared distance on the fiber
            ∇Q = skewlog(Q)[p+1:end, p+1:end]
            ∇V = skewlog(V)
            ∇iter = 0
            while norm(∇Q) + norm(∇V) > ∇tol && ∇iter < ∇max
                tempQ1 .= Q[:, 1:p]
                tempQ2 .= Q[:, p+1:end]
                tempV  .= V
                ∇Q = skewlog(Q)[p+1:end, p+1:end]
                ∇V = skewlog(V)
                G2 =  exp(-2α * ∇Q)
                G1 =  exp(-α * ∇V)
                mul!(Q[:, 1:p], tempQ1, G1, 1, 0)         #Rotate first columns of Q
                mul!(V, tempV, G1, 1, 0)                  #Rotate V   
                mul!(Q[:, p+1:end], tempQ2, G2, 1, 0)     #Rotate last columns of Q  
                ∇iter += 1
            end
            Ω .=  skewlog(Q)
            Ψ .=  skewlog(V)                 
            d = sqrt(β * norm(Ω[1:p, 1:p] - Ψ)^2 + norm(Ω[p+1:end, 1:p])^2)   
            if d < ρ
                print("ρ is greater than the injectivity radius \n")
                return true , iter + 1
            end
        end
        iter +=1
    end
    print("Trial budget exceeded. Either ρ is smaller than the injectivity radius; 
    or it is larger and the algorithm has not been able to find a closer cut time, 
    a likely event when ρ is only slightly larger than the injectivity radius. \n")
    return false , iter
end




