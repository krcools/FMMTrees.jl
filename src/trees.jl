struct Box{T}
    begin_idx::Int
    end_idx::Int
    num_children::Int
    data::T
end

struct ChildView{T}
    tree::T
end

Base.start(itr::ChildView) = (0, 2) # progess, relative_index
Base.done(itr::ChildView, state) = (state[1] == itr.tree[1].num_children)
function Base.next(itr::ChildView, state)
    child = itr.tree[state[2]]
    newstate = (state[1] + child.num_children + 1, state[2] + child.num_children + 1)
    return view(itr.tree, state[2]:endof(itr.tree)) , newstate
end

Base.iteratorsize(cv::ChildView) = SizeUnknown()


children(tree) = ChildView(tree)

function depthfirst(f, t, level = 1)
    f(t, level)
    for c in children(t)
        depthfirst(f, c, level+1)
    end
end

const Tree = AbstractVector{T} where T<:Box
const BlockTree = Tuple{Tree,Tree}
children(b::Tuple{Tree,Tree}) = IterTools.product(children(b[1]), children(b[2]))


function admissable_partition(bt, adm)
    p = Vector{BlockTree}()
    admissable_partition!(bt, adm, p)
    p
end

function admissable_partition!(blocktree, adm, partition)
    adm(blocktree) && (push!(partition, blocktree); return)
    for c in children(blocktree)
        admissable_partition!(c, adm, partition)
    end
end

# function compress(block, adm, partition)
#
#     if adm(blocktree)
#         return lowrank(block)
#     end
#
#     compressed = []
#     for c in children(block)
#         compressed = compressed ∪ compress(c, adm, partition)
#     end
#
#     if length(compressed) == length(children(block))
#         recompressed = recompress(compressed)
#         if size(recompressed) < size(compressed)
#             compressed = recompressed
#         end
#     end
# end


function h1compress(op, tfs, bfs)

    p = positions(tfs)
    q = positions(bfs)

    p2, tp, permp = clustertree(p)
    q2, tq, permq = clustertree(q)
    μ = (τ,σ) -> assemble(op, subset(tfs,τ), subset(bfs,σ))

    η = 1.5
    nmin = 100
    function adm4(b)
        I = b[1][1].begin_idx : b[1][1].end_idx-1
        J = b[2][1].begin_idx : b[2][1].end_idx-1
        length(I) < nmin && return true
        length(J) < nmin && return true
        ll1, ur1 = boundingbox(p2[I]);
        ll2, ur2 = boundingbox(q2[J]);
        c1 = (ll1+ur1)/2;
        c2 = (ll2+ur2)/2;
        diam1 = 2norm(ur1-c1)
        diam2 = 2norm(ur2-c2)
        dist12 = norm(c2-c1) - 0.5(diam1 + diam2)
        return dist12 >= η*min(diam1, diam2)
    end

    block_tree = (tp, tq)
    T = scalartype(op, tfs, bfs)
    hmatrix = Vector{LowRankBlock{T}}()
    h1compress!(block_tree, adm4, hmatrix, μ, permp, permq)

    M = numfunctions(tfs); I = collect(1:M)
    N = numfunctions(bfs); J = collect(1:N)
    return HMatrix(hmatrix, I, J)
end

function h1compress!(block, adm3, hmatrix, μ, permp, permq)

    τ = permp[block[1][1].begin_idx : block[1][1].end_idx-1]
    σ = permq[block[2][1].begin_idx : block[2][1].end_idx-1]

    if adm3(block)
        lrb = aca2(μ,τ,σ)
        #println("Rank: ", size(lrb.matrix.A,2))
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

    #return false

    @assert nc == 4
    if all_leafs
        lrb = aca2(μ,τ,σ)
        sz2 = storedentries(lrb)
        if sz2 <= sz1
            for i in 1:nc; pop!(hmatrix); end
            push!(hmatrix, lrb)
            return true
            println("recompression: $sz1 → $sz2")
        else
            println("recompression failed")
            return false
        end
    end

    return false
end
