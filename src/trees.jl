export print_tree

struct ChildView{T}
    tree::T
end

Base.iteratorsize(cv::ChildView) = SizeUnknown()


children(tree) = ChildView(tree)

function depthfirst(f, t, level = 1)
    f(t, level)
    for c in children(t)
        depthfirst(f, c, level+1)
    end
end

function print_tree(tree)
    depthfirst(tree) do node,level
        print("-"^(level-1))
        print(node[1].data)
        println()
    end
end
