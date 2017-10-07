using IterTools

export clustertree, children, depthfirst, admissable_partition, boundingbox

struct Box
    begin_idx::Int
    end_idx::Int
    sector::Int
end

range(box::Box) = box.begin_idx:box.end_idx-1

function clustertree(points)
    point_list = list(points)

    root_sector = -1 # The root has no well-defined sector
    root_box = Box(1, length(points)+1, root_sector)
    root_node = TreeNode(0, root_box)
    box_list = list([root_node])

    root = start(box_list)
    a = clustertree_impl!(point_list, box_list, root)

    #box_list.data[box_list.nodes[root].value] = TreeNode(a, Box(1, length(points)+1, -1))
    #box_list.data[box_list.node[root].value] = TreeNode(a, root_box)
    box_list[root] = TreeNode(a, root_box)

    # retrieve the permutation linking the original ordering of points
    # and that dictated by its geometric configuration.
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
    end

    return collect(point_list), collect(box_list), perm
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

"""
    clustertree_impl!(points, boxes, boxit)

Insert `points` in the tree backed by `boxes` at the leaf node `boxit`. The function
returns the number of child nodes in the constructed subtree (not counting the root).
"""
function clustertree_impl!(points, boxes, boxit)

    length(points) == 1 && return 0
    box = boxes[boxit]
    @assert length(points)  == box.data.end_idx - box.data.begin_idx

    ll, ur = boundingbox(points)
    diag = ur - ll
    cntr = (ur + ll) / 2

    # divide the pointset along the longest axis of the bb
    i = indmax(diag)
    n1 = count(p -> (p[i] < cntr[i]), points)
    n2 = length(points) - n1

    α, β = box.data.begin_idx, box.data.end_idx

    @assert n1 != 0; @assert n2 != 0
    @assert β - α == n1 + n2

    rbox = TreeNode(0, Box(α+n1, β, 1))
    lbox = TreeNode(0, Box(α, α+n1, 0))

    insert_after!(boxes, rbox, boxit)
    _, rboxit = next(boxes, boxit)

    insert_after!(boxes, lbox, boxit)
    _, lboxit = next(boxes, boxit)

    # move all points in the left box to the front
    head = node = start(points)
    @assert !done(points, node)
    !done(points, node) && (node = advance(points, node))
    while !done(points, node)
        p, newnode = next(points, node)
        if p[i] < cntr[i]
            move_before!(points, node, head)
            #head = node
        end
        node = newnode
    end

    split = start(points)
    for j in 1:n1
        split = advance(points, split)
    end

    @assert (head == split) || (advance(points,head) == split)


    if n1 == 1 && n2 == 1
        a = 2
    elseif n2 == 1
        a1 = clustertree_impl!(sublist(points, start(points), split), boxes, lboxit)
        boxes[lboxit] = TreeNode(a1, Box(α, α+n1, 0))
        a = a1 + 2
    elseif n1 == 1
        a2 = clustertree_impl!(sublist(points, split, done(points)), boxes, rboxit)
        boxes[rboxit] = TreeNode(a2, Box(α+n1, β, 1))
        a = 2 + a2
    else
        a1 = clustertree_impl!(sublist(points, start(points), split), boxes, lboxit)
        boxes[lboxit] = TreeNode(a1, Box(α, α+n1, 0))
        a2 = clustertree_impl!(sublist(points, split, done(points)),  boxes, rboxit)
        boxes[rboxit] = TreeNode(a2, Box(α+n1, β, 1))
        a = 2 + a1 + a2
    end

    a
end
