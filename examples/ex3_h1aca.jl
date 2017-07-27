using FMMTrees, BEAST, CompScienceMeshes

a, h = 1.0, 0.25
Γ1 = meshsphere(a, h)
Γ2 = CompScienceMeshes.translate(Γ1, point(1.8,0,0))

X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)

κ = 1.0; γ = κ*im;
T = BEAST.MWSingleLayer3D(0.0im, 1.0im, 0.0im)

μ = (τ,σ) -> assemble(T, subset(X1,τ), subset(X2,σ))
τ = collect(1:numfunctions(X1))
σ = collect(1:numfunctions(X2))
A = FMMTrees.assemble_aca(T,X1,X2);

Z = assemble(T,X1,X2)
sz = sum(length(a.a)*length(a.row_idcs)+length(a.b)*length(a.col_idcs) for a in A)
println("Compression: ", sz/numfunctions(X1)^2)

Q = FMMTrees.reconstruct(A, numfunctions(X1), numfunctions(X2))
@show norm(Q-Z)
@show norm(Q)
@show norm(Z)
