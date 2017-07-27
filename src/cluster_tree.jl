using IterTools

export clustertree, children, depthfirst, admissable_partition, boundingbox

struct Box{T}
    begin_idx::Int
    end_idx::Int
    num_children::Int
    data::T
end

function clustertree(points)
    point_list = list(points)
    box_list = list([Box{Int}(1, length(points)+1, 0, -1)])
    root = start(box_list)
    a = clustertree_impl!(point_list, box_list, root)
    box_list.data[box_list.nodes[root].value] = Box{Int}(1, length(points)+1, a, -1)

    np = length(point_list)
    perm = zeros(Int, length(point_list))
    points_out = similar(points)

    i = 1
    s = start(point_list)
    while !done(point_list,s)
        @assert 1 <= s-1 <= np
        perm[i] = s-1
        i += 1
        _, s = next(point_list, s)
        #iperm[s] = i
    end

    return collect(point_list), collect(box_list), perm #, iperm
end

function boundingbox(points)
    @assert !isempty(points)
    ll, ur = first(points), first(points)
    for p in points
        ll = min.(ll,p)
        ur = max.(ur,p)
    end
    return ll, ur
end

function clustertree_impl!(points, boxes, boxit)

    length(points) == 1 && return 0
    box = peek(boxes, boxit)
    @assert length(points)  == box.end_idx - box.begin_idx

    ll, ur = boundingbox(points)
    diag = ur - ll
    cntr = (ur + ll) / 2

    # divide the pointset along the longest axis of the bb
    i = indmax(diag)
    n1 = count(p -> (p[i] < cntr[i]), points)
    n2 = length(points) - n1

    α, β = box.begin_idx, box.end_idx
    @assert n1 != 0; @assert n2 != 0
    @assert β - α == n1 + n2
    rbox = Box(α+n1, β, 0, 1); insert_after!(boxes, rbox, boxit); _, rboxit = next(boxes, boxit)
    lbox = Box(α, α+n1, 0, 0); insert_after!(boxes, lbox, boxit); _, lboxit = next(boxes, boxit)

    # move all points in the left node to the front
    head = node = start(points)
    !done(points, node) && ((_, node) = next(points, node))
    while !done(points, node)
        p, newnode = next(points, node)
        if p[i] < cntr[i]
            move_before!(points, node, head)
            head = node
        end
        node = newnode
    end

    split = start(points); for i in 1:n1; _, split = next(points, split); end

    if n1 == 1 && n2 == 1
        a = 2
    elseif n2 == 1
        a1 = clustertree_impl!(sublist(points, start(points), split), boxes, lboxit)
        boxes.data[boxes.nodes[lboxit].value] =  Box(α, α+n1, a1, 0)
        a = a1 + 2
    elseif n1 == 1
        a2 = clustertree_impl!(sublist(points, split, done(points)), boxes, rboxit)
        boxes.data[boxes.nodes[rboxit].value] =  Box(α+n1, β, a2, 1)
        a = 2 + a2
    else
        a1 = clustertree_impl!(sublist(points, start(points), split), boxes, lboxit)
        boxes.data[boxes.nodes[lboxit].value] =  Box(α, α+n1, a1, 0)
        a2 = clustertree_impl!(sublist(points, split, done(points)),  boxes, rboxit)
        boxes.data[boxes.nodes[rboxit].value] =  Box(α+n1, β, a2, 1)
        a = 2 + a1 + a2
    end

    a
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


children(tree) = ChildView(tree)

function depthfirst(f, t, level = 1)
    print("  "^(level-1))
    f(t, level)
    for c in children(t)
        depthfirst(f, c, level+1)
    end
end

const Tree = AbstractVector{T} where T<:Box
const BlockTree = Tuple{Tree,Tree}
children(b::Tuple{Tree,Tree}) = IterTools.product(children(b[1]), children(b[2]))


admissable_partition(bt, adm) = (p = Vector{BlockTree}(); admissable_partition!(bt, adm, p); p)
function admissable_partition!(blocktree, adm, partition)
    adm(blocktree) && (push!(partition, blocktree); return)
    for c in children(blocktree)
        admissable_partition!(c, adm, partition)
    end
end
