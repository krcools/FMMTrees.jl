using StaticArrays
using FMMTrees
using VectorBackedLists

T = Float64
P = SVector{3,T}

num_points = 300
points = rand(P,num_points)

point_list = list(points)

function buildtree(points)

    point_list = list(points)

end


function buildtree_impl(points, box, boxptr)

    world_dim = 3
    num_children = 2^world_dim

    chd_boxes [FMMTrees.PointerBasedTrees.Node(
        (0,0), # placeholder data
        0, # num_children
        boxptr, # parent box
    )

    for p in points
        s = findsector(p, box)
        push!(chil)


end
