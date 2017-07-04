export list, insert_after!, move_before!, prev, peek, sublist

struct Node
    value::Int
    next::Int
    prev::Int
    idx::Int
end

struct List{T,S<:AbstractVector{T}}
    data::S
    nodes::Vector{Node}
    head::Int
    tail::Int
end

#value(list, node) = list.data[node.value]

Base.eltype{T,S}(::Type{List{T,S}}) = T
Base.start(list::List) = list.nodes[1].next
Base.next(list::List, state) = (list.data[list.nodes[state].value], list.nodes[state].next)
Base.done(list::List, state) = (list.nodes[state].value == 0)
Base.done(list::List) = list.tail
Base.length(list::List) = length(list.data)

Base.setindex!(list::List, state, v) = (list.data[list.nodes[state].value] = v)
Base.getindex(list::List, state) = list.data[list.nodes[state].value]

peek(list, state) = list.data[list.nodes[state].value]

# sublist iteration
struct SubList{L<:List}
    parent::L
    head::Int
    done::Int
end

sublist(itr, b, e) = SubList(itr, itr.nodes[itr.nodes[b].prev].idx, e)
sublist(ls::SubList, b, e) = sublist(ls.parent, b, e)

Base.start(sl::SubList) = next(sl.parent, sl.head)[2]
Base.next(sl::SubList, s::Int) = next(sl.parent, s)
Base.done(sl::SubList, s::Int) = (s == sl.done)
Base.done(sl::SubList) = sl.done
Base.length(sl::SubList) = (n = 0; for x in sl; n += 1; end; n)

prev(list, state) = (list.data[list.nodes[state].value], list.nodes[state].prev)

function list(data)
    n = length(data)
    nodes = Vector{Node}(n+2)
    nodes[1] = Node(0,2,0,1)
    for i in 2:n+1; nodes[i] = Node(i-1, i+1, i-1, i); end
    nodes[end] = Node(0,0,n+1,n+2)
    List{eltype(data), typeof(data)}(data, nodes, 1, n+2)
end

# i: moving node
# t: new successor
function move_before!(list, I, T)

    @assert I != T
    nodes = list.nodes

    # step 1: remove n
    _, P = prev(list,I);
    _, N = next(list,I);
    p = nodes[P]
    n = nodes[N]
    @assert P == p.idx
    @assert N == n.idx

    nodes[P] = Node(p.value, n.idx, p.prev, p.idx)
    nodes[N] = Node(n.value, n.next, p.idx, n.idx)

    # step 2: reintroduce n
    _, Q = prev(list, T)
    i = nodes[I]
    t = nodes[T]
    q = nodes[Q]
    @assert Q == q.idx
    @assert T == t.idx
    @assert I == i.idx

    nodes[Q] = Node(q.value, i.idx, q.prev, q.idx)
    nodes[T] = Node(t.value, t.next, i.idx, t.idx)
    nodes[I] = Node(i.value, t.idx, q.idx, i.idx)
    nothing
end

move_before!(ls::SubList, item, target) = move_before!(ls.parent, item, target)

function insert_after!(list::List, v, T)

    data = list.data
    nodes = list.nodes

    push!(data, v)

    _, N = next(list, T)
    t = nodes[T]
    n = nodes[N]

    I = length(nodes)+1
    #push!(nodes, Node(length(data), t.next, N, I))
    push!(nodes, Node(length(data), N, T, I))

    nodes[T] = Node(t.value, I, t.prev ,T)
    nodes[N] = Node(n.value, n.next, I, N)

    nothing
end

insert_after!(ls::SubList, value, target) = insert_after!(ls.parent, value, target)
