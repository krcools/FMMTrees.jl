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

function compress(block, adm, partition)

    if adm(blocktree)
        return lowrank(block)
    end

    compressed = []
    for c in children(block)
        compressed = compressed ∪ compress(c, adm, partition)
    end

    if length(compressed) == length(children(block))
        recompressed = recompress(compressed)
        if size(recompressed) < size(compressed)
            compressed = recompressed
        end
    end
end
