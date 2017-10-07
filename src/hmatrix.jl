export reconstruct, storedentries

struct LowRankMatrix{T}
    A::Matrix{T}
    B::Matrix{T}
end

Base.:*(M::LowRankMatrix, x::Vector) = (M.A*(M.B*x))

struct LowRankBlock{T}
    matrix::LowRankMatrix{T}
    τ::Vector{Int}
    σ::Vector{Int}
end

struct HMatrix{T}
    blocks::Vector{LowRankBlock{T}}
    I::Vector{Int}
    J::Vector{Int}
end

Base.eltype(A::HMatrix{T}) where {T} = T
Base.size(A::HMatrix) = (length(A.I), length(A.J))

function Base.:*(A::HMatrix, x)
    T = promote_type(eltype(A), eltype(x))
    y = zeros(T, length(A.I))
    for block in A.blocks
        y[block.τ] .+= block.matrix * x[block.σ]
    end
    return y
end



storedentries(lrm::LowRankMatrix) = length(lrm.A) + length(lrm.B)
storedentries(lrb::LowRankBlock) = storedentries(lrb.matrix)

function storedentries(A::HMatrix)
    sz = 0
    for block in A.blocks
        sz += storedentries(block)
    end
    return sz
end

function reconstruct(A::HMatrix)
    M = zeros(eltype(A), size(A)...)
    for block in A.blocks
        @assert norm(M[block.τ,block.σ]) == 0
        M[block.τ,block.σ] = block.matrix.A * block.matrix.B
    end
    return M
end


function admissability_condition(strength, nmin, test_points, trial_points)

    η = strength

    p2 = test_points
    q2 = trial_points

    function adm(b)

        τ = range(data(testcluster(b)))
        σ = range(data(trialcluster(b)))

        length(τ) < nmin && return true
        length(σ) < nmin && return true

        ll1, ur1 = boundingbox(p2[τ]);
        ll2, ur2 = boundingbox(q2[σ]);

        c1 = (ll1+ur1)/2;
        c2 = (ll2+ur2)/2;

        diam1 = 2norm(ur1-c1)
        diam2 = 2norm(ur2-c2)

        dist12 = norm(c2-c1) - 0.5(diam1 + diam2)
        return dist12 >= η*min(diam1, diam2)
    end
end


function h1compress(op, tfs, bfs)

    T = scalartype(op, tfs, bfs)

    p = positions(tfs)
    q = positions(bfs)

    p2, tp, permp = clustertree(p)
    q2, tq, permq = clustertree(q)
    μ = (τ,σ) -> assemble(op, subset(tfs,τ), subset(bfs,σ))

    assembler = blockassembler(op, tfs, bfs)
    function μ2(τ,σ)
        Z = zeros(T,length(τ),length(σ))
        assembler(τ,σ,(v,m,n)->(Z[m,n] += v))
        return Z
    end

    η, nmin = 1.5, 100
    adm4 = admissability_condition(η, nmin, p2, q2)

    block_tree = (tp, tq)
    T = scalartype(op, tfs, bfs)
    hmatrix = Vector{LowRankBlock{T}}()
    h1compress!(block_tree, adm4, hmatrix, μ2, permp, permq)

    M = numfunctions(tfs)
    N = numfunctions(bfs)

    I = collect(1:M)
    J = collect(1:N)

    return HMatrix(hmatrix, I, J)
end


function h1compress!(block, adm3, hmatrix, μ, permp, permq)

    τ = permp[range(data(testcluster(block)))]
    σ = permq[range(data(trialcluster(block)))]

    if adm3(block)
        lrb = aca2(μ,τ,σ)
        rank1 = min(length(τ), length(σ))
        rank2 = size(lrb.matrix.A,2)
        #println("Compression: $rank1 ⇢ $rank2")
        push!(hmatrix, lrb);
        return true
    end

    all_leafs = true
    sz1, nc = 0, 0
    for child in children(block)
        nc += 1
        is_leaf = h1compress!(child, adm3, hmatrix, μ, permp, permq)
        if is_leaf
            sz1 += storedentries(last(hmatrix))
        else
            all_leafs = false
        end
    end

    # WARNING: Recompression disabled
    return false

    # @assert nc == 4
    # if all_leafs
    #     lrb = aca2(μ,τ,σ)
    #     sz2 = storedentries(lrb)
    #     if sz2 <= sz1
    #         for i in 1:nc; pop!(hmatrix); end
    #         push!(hmatrix, lrb)
    #         return true
    #         println("recompression: $sz1 → $sz2")
    #     else
    #         println("recompression failed")
    #         return false
    #     end
    # end
    #
    # return false
end
