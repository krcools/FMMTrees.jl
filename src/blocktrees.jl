const BlockTree = Tuple{VectorBackedTree,VectorBackedTree}
testcluster(blocktree::BlockTree) = blocktree[1]
trialcluster(blocktree::BlockTree) = blocktree[2]
children(b::BlockTree) = IterTools.product(children(testcluster(b)), children(trialcluster(b)))


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
