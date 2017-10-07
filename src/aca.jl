using BEAST
using DataStructures

function indmax2(f, itr)
    J = first(itr)
    F = f(J)
    for i in itr
        fi = f(i)
        fi > F && (J = i; F = fi)
    end
    return J, F
end

# μ is a lazy array. After all, the whole point of doing ACA is to
# avoid explicit computation of all entries of the matrix under consideration.
function adaptive_cross_approximation(μ,τ,σ,T=Complex128)

    ϵ = sqrt(eps(real(T)))*1000

    @assert !isempty(τ)
    @assert !isempty(σ)

    I = OrderedSet(τ)
    J = OrderedSet(σ)

    αs = Vector{T}()
    as = Vector{Vector{T}}()
    bs = Vector{Vector{T}}()

    l, q, j = 0, 1, first(J)
    delete!(J,j)
    @assert σ[1] == j

    R1ic = zeros(T, length(σ))
    R1cj = zeros(T, length(τ))

    Ric = zeros(T, length(σ))
    Rcj = zeros(T, length(τ))

    while true
        l += 1

        fill!(Rcj,0)
        for k in eachindex(αs)
            Rcj .+= αs[k]*as[k]*bs[k][q]
        end
        Mcj = vec(μ(τ,[j]))

        isempty(I) && break

        Iv = collect(I) # I is a subset of tau
        P = indexin(Iv,τ)
        p = indmax(abs.(Mcj[P] - Rcj[P]))
        i = Iv[p]
        p = P[p]

        @assert !(norm(Mcj[P] - Rcj[P]) + norm(Mcj[P]) ≈ norm(Mcj[P]))
        @assert !( abs(Mcj[p]-R1cj[p]) + abs(Mcj[p]) ≈ abs(Mcj[p]))
        @assert τ[p] == i
        @assert 1 <= p <= length(τ)

        delete!(I,i)
        fill!(Ric,0)
        for k in eachindex(αs)
            Ric .+= αs[k]*as[k][p]*bs[k]
        end

        Mic = vec(μ([i],σ))

        @assert Rcj[p] ≈ Ric[q]
        @assert Mcj[p] ≈ Mic[q]
        @assert !(norm(Mic-Ric) + norm(Mic) ≈ norm(Mic))
        α = 1/(Mic[q]-Ric[q])
        @assert !isnan(α)
        a = Mcj-Rcj
        b = Mic-Ric

        l > 1 && norm(α)*norm(a)*norm(b) < ϵ && break

        push!(αs, α)
        push!(as, a)
        push!(bs, b)

        isempty(J) && break
        Jv = collect(J)
        Q = indexin(Jv,σ)
        q = indmax(abs.(Mic[Q] - R1ic[Q]))
        j = Jv[q]
        q = Q[q]
        @assert σ[q] == j
        @assert 1 <= q <= length(σ)
        delete!(J,j)

        R1ic = Ric
        R1cj = Rcj
    end

    return αs, as, bs, l
end


function aca2(μ,τ,σ,T=Complex128)
    cs, as, bs, l = adaptive_cross_approximation(μ,τ,σ,T)
    r = length(cs)
    m,n = length(τ), length(σ)
    A, B = zeros(T,m,r), zeros(T,r,n)
    for i in 1:r; A[:,i] .= cs[i]*as[i]; end
    for i in 1:r; B[i,:] .= bs[i]; end
    lrm = LowRankMatrix(A,B)
    #lrm = recompress2(lrm)
    return LowRankBlock(lrm,τ,σ)
end


function recompress2(lrm, ϵ = sqrt(eps(real(eltype(lrm.A))))*1000)

    A, B = lrm.A, lrm.B
    r = size(A,2)

    Q,R = qr(A)
    U,s,V = svd(R*B)
    r′ = findfirst(x -> abs(x) < ϵ ,s)
    r′ == 0 && (r′ = r)

    A = (Q*U)[:,1:r′]
    B = (diagm(s)*V')[1:r′,:]

    #println("Compression: $r → $r′")

    return LowRankMatrix(A,B)
end

function recompress(αs,as,bs,ϵ = sqrt(eps(real(eltype(αs))))*100)

    @assert !isempty(αs)
    T = eltype(αs)

    r = length(αs)
    m = length(as[1])
    n = length(bs[1])

    A = zeros(T,m,r)
    B = zeros(T,r,n)
    for i ∈ 1:r
        A[:,i] = αs[i] * as[i]
        B[i,:] = bs[i]
    end

    Q,R = qr(A)
    U,s,V = svd(R*B)
    @assert size(s) == (r,)
    r′ = findfirst(x -> abs(x) < ϵ ,s)
    r′ == 0 && (r′ = r)

    A = (Q*U)
    B = diagm(s)*V'

    as = [A[:,i] for i in 1:r′]
    bs = [B[i,:] for i in 1:r′]
    αs = ones(T,r′)

    return αs, as, bs, r′
end
