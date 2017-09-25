using FMMTrees, BEAST, CompScienceMeshes

a, h = 1.0, 0.061
Γ1 = meshsphere(a, h)
Γ2 = CompScienceMeshes.translate(Γ1, point(2.10,0,0))

X1 = raviartthomas(Γ1)
X2 = raviartthomas(Γ2)

@show numfunctions(X1)
@show numfunctions(X2)

κ = 1.0; γ = κ*im;
T = Maxwell3D.singlelayer(gamma=γ)

μ = (τ,σ) -> assemble(T, subset(X1,τ), subset(X2,σ))
#A = FMMTrees.assemble_aca(T,X1,X2);
A = FMMTrees.h1compress(T,X1,X2)



sz = storedentries(A)
println("Compression: ", sz/numfunctions(X1)/numfunctions(X2))

# Z = assemble(T,X1,X2)
# Q = reconstruct(A)
# @show norm(Q-Z)
# @show norm(Q)
# @show norm(Z)


# p1 = positions(X1)
# q1, t1, perm1 = clustertree(p1)
# iperm1 = collect(1:numfunctions(X1))
# ipermute!(iperm1, perm1)
#
# p2 = positions(X2)
# q2, t2, perm2 = clustertree(p2)
# iperm2 = collect(1:numfunctions(X2))
# ipermute!(iperm2, perm2)
#
#
# B = zeros(numfunctions(X1),numfunctions(X2))
# for b in A.blocks
#     τ, σ = iperm1[b.τ], iperm2[b.σ]
#     B[τ,σ] .+= log10(length(τ)*length(σ))
# end
