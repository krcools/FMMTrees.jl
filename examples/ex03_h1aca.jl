using FMMTrees, BEAST, CompScienceMeshes

#a, h = 1.0, 0.061
a, h = 1.0, 0.151
Γ1 = meshsphere(a, h)
Γ2 = CompScienceMeshes.translate(Γ1, point(2.10,0,0))

X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)

@show numfunctions(X1)
@show numfunctions(X2)

κ = 1.0; γ = κ*im;
T = Maxwell3D.singlelayer(gamma=γ)

μ = (τ,σ) -> assemble(T, subset(X1,τ), subset(X2,σ))
A = FMMTrees.h1compress(T,X1,X2)

sz = storedentries(A)
println("Compression: ", sz/numfunctions(X1)/numfunctions(X2))

# Z1 = assemble(T,X1,X2)
# Z2 = FMMTrees.reconstruct(A)
