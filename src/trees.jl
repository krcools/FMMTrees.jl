export data

"""
    data(node)

Retrieve the data aka payload associated with the given node.
"""
function data end

# import AbstractTrees: print_tree, children
# export print_tree, children
#
# export insertchild, data
#
# function insertchild end
# function data end
#
#
# function depthfirst(f, t, level = 1)
#     f(t, level)
#     for c in children(t)
#         depthfirst(f, c, level+1)
#     end
# end

# function print_tree(tree)
#     depthfirst(tree) do tree, level
#         print("-"^(level-1))
#         print(data(tree))
#         println()
#     end
# end
