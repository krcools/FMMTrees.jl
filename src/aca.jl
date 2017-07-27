using BEAST
using DataStructures

function indmax2(f, itr)
    I = first(itr)
    F = f(I)
    for i in itr
        fi = f(i)
        fi > F && (I = i; F = fi)
    end
    return I, F
end

# μ is a lazy array. After all, the whole point of doing ACA is to
# avoid explicit computation of all entries of the matrix under consideration.
function adaptive_cross_approximation(μ,τ,σ)

    T = Complex128 # TODO: generalise
    ϵ = sqrt(eps(real(T)))*100

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
        #@show l

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

        norm(α)*norm(a)*norm(b) < ϵ && break

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

    return αs, as, bs
end

# dop: discrete operator (operator, testfunctions, trialfunctions) triple
function assemble_aca(op, tfs, bfs)

    p = positions(tfs)
    q = positions(bfs)

    p, tp, permp = clustertree(p)
    q, tq, permq = clustertree(q)

    η = 2.0
    @show η

    # create closure
    function adm2(b)
        nmin = 20
        #η = 1.5
        #η = 2.0
        #η = 2.5
        I = b[1][1].begin_idx : b[1][1].end_idx-1
        J = b[2][1].begin_idx : b[2][1].end_idx-1
        length(I) < nmin && return true
        length(J) < nmin && return true
        #ll1, ur1 = boundingbox(p[permp[I]]); c1 = (ll1+ur1)/2;
        #ll2, ur2 = boundingbox(q[permq[J]]); c2 = (ll2+ur2)/2;
        ll1, ur1 = boundingbox(p[I]); c1 = (ll1+ur1)/2;
        ll2, ur2 = boundingbox(q[J]); c2 = (ll2+ur2)/2;
        diam1 = norm(ur1-c1)
        diam2 = norm(ur2-c2)
        dist12 = norm(c2-c1)
        # @show diam1
        # @show diam2
        # @show dist12
        # error("basta")
        return dist12 >= η*max(diam1, diam2)
    end

    P = admissable_partition((tp,tq), adm2)
    @show length(P)
    μ = (τ,σ) -> assemble(op, subset(tfs,τ), subset(bfs,σ))

    T = scalartype(op)
    A = Vector{LowRankBlock{T}}()
    for p in P
        τ = p[1][1].begin_idx : p[1][1].end_idx-1
        σ = p[2][1].begin_idx : p[2][1].end_idx-1
        α, a, b = adaptive_cross_approximation(μ, permp[τ], permq[σ])
        push!(A, LowRankBlock(α,a,b,permp[τ],permq[σ]))
    end

    return A
end

struct LowRankBlock{T} #<: AbstractMatrix{T}
    α::Vector{T}
    a::Vector{Vector{T}}
    b::Vector{Vector{T}}
    row_idcs
    col_idcs
end

function reconstruct(aca::Vector{M}, m, n) where M <: LowRankBlock
    T = eltype(aca[1].α)
    A = zeros(T,m,n)
    for b in aca
        R = sum(c*a*b.' for (c,a,b) in zip(b.α,b.a,b.b))
        A[b.row_idcs, b.col_idcs] .+= R
    end
    A
end
export reconstruct
