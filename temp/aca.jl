using BEAST, CompScienceMeshes

m1 = meshcuboid(1.0, 1.0, 1.0, 0.05)
m2 = meshcuboid(1.0, 1.0, 1.0, 0.05)
CompScienceMeshes.translate!(m2, point(20,0,0))

X1 = raviartthomas(m1)
X2 = raviartthomas(m2)
r = numfunctions(X1) / 2

ts = BEAST.MWSingleLayer3D(0.0, 1.0, 0.0)
M = assemble(ts,X1,X2)

## create the low rank approximation

T = Float64
αs = T[]
as = Vector{T}[]
bs = Vector{T}[]

m, n = size(M)
I = Set(1:m)
J = Set(1:n)

r = 30
l = 0
R = zeros(T,m,n)
R1 = zeros(T,m,n)
R2 = zeros(T,m,n)
j = first(J)
delete!(J,j)
@assert j ∉ J
while true
    l += 1
    @show l
    l > r && break
    @show maximum(abs(M[k,j] - R[k,j]) for k in 1:m)
    i = indmax(abs(M[k,j] - R[k,j]) for k in 1:m)
    delete!(I,i)
    @assert i ∉ I
    α, a, b = 1/(M[i,j]-R[i,j]), M[:,j]-R[:,j], M[i,:]-R[i,:]
    @assert M[i,j]-R[i,j] != 0

    push!(αs, α)
    push!(as, a)
    push!(bs, b)

    @show maximum(abs(M[i,k] - R1[i,k]) for k in J)
    j = indmax(abs(M[i,k] - R1[i,k]) for k in J)
    delete!(J,j)
    @assert j ∉ J

    R2 = R1
    R1 = R
    R += α*a*b'
end

function matvec2(αs, as, bs, x, y)
    for i in eachindex(as)
        α, a, b = αs[i], as[i], bs[i]
        y .+= α * dot(b,x) * a
    end
end
