module BlockTrees

import FMMTrees

struct BlockTree{T}
    test_cluster::T
    trial_cluster::T
end

testcluster(blocktree::BlockTree) = blocktree.test_cluster
trialcluster(blocktree::BlockTree) = blocktree.trial_cluster

FMMTrees.root(t::BlockTree) = (FMMTrees.root(testcluster(t)), FMMTrees.root(trialcluster(t)))
function FMMTrees.data(t::BlockTree, n)
    (
        FMMTrees.data(testcluster(t), n[1]),
        FMMTrees.data(trialcluster(t), n[2]))
end

function FMMTrees.children(b::BlockTree, node)
    test_chds = FMMTrees.children(testcluster(b), node[1])
    trial_chds = FMMTrees.children(testcluster(b), node[2])
    ((ch[1],ch[2]) for ch in Iterators.product(test_chds, trial_chds))
end

function FMMTrees.haschildren(b::BlockTree, node)
    !FMMTrees.haschildren(testcluster(b), node[1]) && return false
    !FMMTrees.haschildren(trialcluster(b), node[2]) && return false
    return true
end

end # module BlockTrees
