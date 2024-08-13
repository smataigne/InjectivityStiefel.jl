using LinearAlgebra, SkewLinearAlgebra
#=
a = 1im*sqrt(6)/6
b = sqrt(2)/2
c = sqrt(3)/6 + 1im/2
d = sqrt(3)/3

V = [-a a -d -d;
    -a a  c   c';
    -a a  c'  c;
    b  b  0   0]
D1 = Diagonal(1im .*[-π; π; 2π; -2π])
D2 = Diagonal(1im .*[-π; π; 0; 0])
M1 = real.(V * D1* V')
M2 = real.(V * D2* V')
#display(exp(0.3*M1))
#display(exp(0.3*M2))
A = π * sqrt(3)/3 *[0 1 -1;-1 0 1;1 -1 0]
d = sqrt(norm(M1[1:3,1:3] + A)^2 + norm(M1[4,1:3])^2)
#display(d)
#display(π*sqrt(3))
=#
β = 2/3+0.05
A = [0 π/(1-2β); -π/(1-2β) 0]
B = randn(1, 2)
B ./= norm(B)
B .*= 2π * sqrt(1-β^2/(1-2β)^2)
M = [2β*A -B';B 0]
display(exp(M)[:,1:2]*exp((1-2β)*A))
#display(norm((1-2β)*A)*sqrt(2))
#display(π * sqrt(2β^2/(1-2β)^2 + 4*(1-β^2/(1-2β)^2)))

#=
D = randn(2, 2)
D = D - D'
I2 = Matrix(1.0I, 2, 2)
O2 = zeros(2,2)
Ω = skewhermitian!(copy([O2 -I2; I2 O2]))
T = skewhermitian!(copy([D O2;O2 O2]))
Ψ = skewhermitian!(copy(0.5 * D))
ε = 1e-6
t = 0.5
display(exp(t*Ω) * T)
display(T)
f(ε) =  (exp(t * Ω + ε * T)[:,1:2] * exp(-ε * Ψ) - exp(t * Ω )[:,1:2])/ε
display(round.(f(ε), digits = 4))
display( (exp(t * Ω + ε * T)[:,1:2] - exp(t * Ω)[:,1:2])/ε  - exp(t * Ω )[:,1:2] * Ψ)
display([0.5sin(t)/t* skewhermitian!(copy(D));0.5*sin(t) * skewhermitian!(copy(D))])
=#