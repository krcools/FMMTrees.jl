using FMMTrees, BEAST, CompScienceMeshes

a, h = 1.0, 0.11
Γ1 = meshsphere(a, h)
Γ2 = CompScienceMeshes.translate(Γ1, point(1.95,0,0))

X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)

κ = 1.0; γ = κ*im;
T = Maxwell3D.singlelayer(gamma=γ)

μ = (τ,σ) -> assemble(T, subset(X1,τ), subset(X2,σ))
#A = FMMTrees.assemble_aca(T,X1,X2);
A = FMMTrees.h1compress(T,X1,X2)

sz = storedentries(A)
println("Compression: ", sz/numfunctions(X1)/numfunctions(X2))

Z = assemble(T,X1,X2)
Q = reconstruct(A)
@show norm(Q-Z)
@show norm(Q)
@show norm(Z)
