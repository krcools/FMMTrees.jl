using FMMTrees, BEAST, CompScienceMeshes

a, h = 1.0, 0.07
Γ1 = Γ2 = meshsphere(a, h)

X1 = X2 = raviartthomas(Γ1)

T = Maxwell3D.singlelayer(gamma=1.0im, alpha=1.0, beta=1.0)
#μ = (τ,σ) -> assemble(T, subset(X1,τ), subset(X2,σ))

A = FMMTrees.assemble_aca(T,X1,X2);

#Z = assemble(T,X1,X2)
sz = storedentries(A)
println("Compression: ", sz/numfunctions(X1)/numfunctions(X2))

# Q = reconstruct(A)
# @show norm(Q-Z)
# @show norm(Q)
# @show norm(Z)
