export list, insert_after!, move_before!, prev, peek, sublist

struct Node
    value::Int  # index into value of the accompanying value
    next::Int   # index into nodes of the next node
    prev::Int   # index into nodes of the previous node
    idx::Int    # index into nodes of the current node
end

struct List{T,S<:AbstractVector{T}}
    data::S
    nodes::Vector{Node}
    head::Int
    tail::Int
end

Base.eltype(::Type{List{T,S}}) where {S,T} = T
Base.start(list::List) = list.nodes[1].next
Base.next(list::List, state) = (list.data[list.nodes[state].value], list.nodes[state].next)
Base.done(list::List, state) = (list.nodes[state].value == 0)

"""
    done(iterable) -> state

Produces an iterator state for which `done(iterable,state) == true`. Cf. to the
c++ end() api.
"""
Base.done(list::List) = list.tail
Base.length(list::List) = length(list.data)

Base.setindex!(list::List, v, state) = (list.data[list.nodes[state].value] = v)
Base.getindex(list::List, state) = list.data[list.nodes[state].value]

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

"""
    prev(list, state) -> item, prevstate

Returns the current item from `list` and sets the state to point to the previous
entry. It hold that

```
_, p = prev(list, s)
_, n = next(list, p)
n == s
```
"""
prev(list, state) = (list.data[list.nodes[state].value], list.nodes[state].prev)

"""
Create a list from an indexable container. The list provided a view on the container,
so any mutations realised through calling the list API will be reflected in the
state of the underlying container.
"""
function list(data)
    n = length(data)
    nodes = Vector{Node}(n+2)
    nodes[1] = Node(0,2,0,1)
    for i in 2:n+1; nodes[i] = Node(i-1, i+1, i-1, i); end
    nodes[end] = Node(0,0,n+1,n+2)
    List{eltype(data), typeof(data)}(data, nodes, 1, n+2)
end

"""
    move_before(list, item, dest)

Move the value pointed to by iterator `item` in fron of iterator `state`.
"""
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

"""
    insert_after!(list, value, dest)

Insert `value` in `list` after the value pointed to by iterator `dest`.
"""
function insert_after!(list::List, v, T)

    data = list.data
    nodes = list.nodes

    push!(data, v)

    _, N = next(list, T)
    t = nodes[T]
    n = nodes[N]

    I = length(nodes)+1
    push!(nodes, Node(length(data), N, T, I))

    nodes[T] = Node(t.value, I, t.prev ,T)
    nodes[N] = Node(n.value, n.next, I, N)
    nothing
end

insert_after!(ls::SubList, value, target) = insert_after!(ls.parent, value, target)
