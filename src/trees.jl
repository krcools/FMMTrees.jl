export print_tree

export root, children, insertchild

function root end
function children end
function insertchild end

struct ChildView{T,N}
    tree::T
    node::N
end

Base.iteratorsize(cv::ChildView) = SizeUnknown()


children(tree,node) = ChildView(tree,node)

function depthfirst(f, t, n, level = 1)
    f(t, n, level)
    for c in children(t, n)
        depthfirst(f, c, level+1)
    end
end

function print_tree(tree)

    depthfirst(tree, root(tree)) do node,level
        print("-"^(level-1))
        print(node[1].data)
        println()
    end
end
