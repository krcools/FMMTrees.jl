export print_tree

struct TreeNode{T}
    num_children::Int
    # begin_idx::Int
    # end_idx::Int
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

const Tree = AbstractVector{T} where T<:TreeNode
data(tree::Tree) = tree[1].data

const BlockTree = Tuple{Tree,Tree}
testcluster(blocktree::BlockTree) = blocktree[1]
trialcluster(blocktree::BlockTree) = blocktree[2]
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


function print_tree(tree)
    depthfirst(tree) do node,level
        print("-"^(level-1))
        print(node[1].begin_idx,",",node[1].end_idx)
        println()
    end
end
