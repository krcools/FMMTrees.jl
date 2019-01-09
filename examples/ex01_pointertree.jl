using FMMTrees
using StaticArrays

num_points

T = Float64
P = SVector{3,T}
points = rand(P, num_points)

# Build a pointer bases tree on top of these points
function buildtree(points; min_box_size=0.1)

    bb = boundingbox(pts)

end


function subdivide(points, box, min_box_size)

    num_child_boxes = 4

    diameter(box) <= min_box_size && return

    C = typeof(points)
    points_partition = [C() for _ in 1:num_child_boxes]
    for point in points
        sector = determine_sector(point, box)
        push!(points_partition[sector], point)
    end

    for p in points_partition
        push!(box.children, Box(points))
    end

end
