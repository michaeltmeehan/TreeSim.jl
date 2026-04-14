import Base: parent

"""
    NodeKind

Structural role of a node in the canonical `Tree` representation.

- `Root`: the unique node with `parent == 0`; it must have at least one child.
- `Binary`: an internal non-root node with exactly two children.
- `UnsampledUnary`: an unsampled internal non-root node with exactly one child.
- `SampledUnary`: a sampled internal non-root node with exactly one child.
- `SampledLeaf`: a sampled non-root tip with no children.
"""
@enum NodeKind::UInt8 begin
    Root = 0
    Binary = 1
    UnsampledUnary = 2
    SampledUnary = 3
    SampledLeaf = 4
end

"""
    Tree

Canonical time-ordered tree storage.

Nodes are identified by their vector index. The integer value `0` is the only
null reference and denotes the absence of a parent or child. Every non-zero
entry in `left`, `right`, and `parent` must point to another node index in the
same `Tree`.

Canonical trees are ordered by nondecreasing node time in storage, and every
edge must point from an earlier index to a later index with strictly increasing
time. A tree that passes `validate_tree` has one root, is reachable from that
root, has reciprocal parent/child pointers, and has node degrees compatible
with `kind`.
"""
struct Tree
    time::Vector{Float64}
    left::Vector{Int}
    right::Vector{Int}
    parent::Vector{Int}
    kind::Vector{NodeKind}
    host::Vector{Int}
    label::Vector{Int}
end

"""
    Node

Immutable row view returned by `tree[i]`. The `id` field is the node's index in
the backing `Tree`; other fields mirror the corresponding vector entries.
"""
struct Node
    id::Int
    time::Float64
    left::Int
    right::Int
    parent::Int
    kind::NodeKind
    host::Int
    label::Int
end

Tree() = Tree(Float64[], Int[], Int[], Int[], NodeKind[], Int[], Int[])

Base.length(tree::Tree) = length(tree.time)
Base.size(tree::Tree) = (length(tree),)
Base.eachindex(tree::Tree) = eachindex(tree.time)

function Base.getindex(tree::Tree, i::Int)
    _check_node(tree, i)
    return Node(i,
                tree.time[i],
                tree.left[i],
                tree.right[i],
                tree.parent[i],
                tree.kind[i],
                tree.host[i],
                tree.label[i])
end

function Base.iterate(tree::Tree, i::Int=1)
    i > length(tree) && return nothing
    return (tree[i], i + 1)
end

Base.pairs(tree::Tree) = zip(eachindex(tree), tree)

function _check_node(tree::Tree, i::Int)
    1 <= i <= length(tree) || throw(BoundsError(tree, i))
    return nothing
end

function children(tree::Tree, i::Int)
    _check_node(tree, i)
    out = Int[]
    tree.left[i] != 0 && push!(out, tree.left[i])
    tree.right[i] != 0 && push!(out, tree.right[i])
    return out
end

function parent(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.parent[i]
end

function isleaf(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.left[i] == 0 && tree.right[i] == 0
end

isinternal(tree::Tree, i::Int) = !isleaf(tree, i)

function isbinary(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.left[i] != 0 && tree.right[i] != 0
end

function isunary(tree::Tree, i::Int)
    _check_node(tree, i)
    return (tree.left[i] != 0) ⊻ (tree.right[i] != 0)
end

function leaves(tree::Tree)
    return [i for i in eachindex(tree) if isleaf(tree, i)]
end

function ancestors(tree::Tree, i::Int)
    _check_node(tree, i)
    out = Int[]
    seen = Set{Int}()
    current = tree.parent[i]
    while current != 0
        _check_node(tree, current)
        current in seen && error("Cycle detected in parent chain from node $i.")
        push!(seen, current)
        push!(out, current)
        current = tree.parent[current]
    end
    return out
end

function descendants(tree::Tree, i::Int)
    _check_node(tree, i)
    out = Int[]
    stack = reverse(children(tree, i))
    seen = Set{Int}()
    while !isempty(stack)
        current = pop!(stack)
        _check_node(tree, current)
        current in seen && error("Cycle detected below node $i.")
        push!(seen, current)
        push!(out, current)
        append!(stack, reverse(children(tree, current)))
    end
    return out
end

function branch_length(tree::Tree, i::Int)
    _check_node(tree, i)
    p = tree.parent[i]
    p == 0 && error("Root node $i has no branch length.")
    _check_node(tree, p)
    return tree.time[i] - tree.time[p]
end

"""
    validate_tree(tree::Tree; require_single_root::Bool=true, require_reachable::Bool=true)

Validate the canonical structural contract for `Tree`.

By default, validation is intentionally strict: non-empty trees must have
exactly one root and all nodes must be reachable from it. Storage order is also
part of the contract: `time` must be nondecreasing by node index, while each
parent-child edge must advance to a greater index and a strictly later time.
Set `require_single_root=false` or `require_reachable=false` only when
inspecting provisional, non-canonical structures.
"""
function validate_tree(tree::Tree; require_single_root::Bool=true, require_reachable::Bool=true)
    n = length(tree.time)
    _validate_lengths(tree, n)
    _validate_indices(tree, n)
    roots = _root_nodes(tree)
    _validate_roots(tree, roots; require_single_root)
    _validate_unique_children(tree)
    _validate_parent_child_consistency(tree)
    _validate_degrees(tree)
    _validate_time_ordering(tree)
    _validate_acyclic(tree)
    require_reachable && _validate_reachable(tree, roots)
    return true
end

function _validate_lengths(tree::Tree, n::Int)
    if !(length(tree.left) == n &&
         length(tree.right) == n &&
         length(tree.parent) == n &&
         length(tree.kind) == n &&
         length(tree.host) == n &&
         length(tree.label) == n)
        error("Tree vectors must have equal lengths.")
    end
    return nothing
end

function _validate_indices(tree::Tree, n::Int)
    for i in 1:n
        for (name, j) in (("left", tree.left[i]), ("right", tree.right[i]), ("parent", tree.parent[i]))
            if j != 0 && !(1 <= j <= n)
                error("Invalid $name index $j at node $i.")
            end
            j == i && error("Node $i cannot reference itself as $name.")
        end
    end
    return nothing
end

function _validate_parent_child_consistency(tree::Tree)
    for i in eachindex(tree)
        l = tree.left[i]
        r = tree.right[i]
        l != 0 && tree.parent[l] != i && error("Parent-child mismatch: node $i left child $l has parent $(tree.parent[l]).")
        r != 0 && tree.parent[r] != i && error("Parent-child mismatch: node $i right child $r has parent $(tree.parent[r]).")
    end

    for i in eachindex(tree)
        p = tree.parent[i]
        if p != 0 && tree.left[p] != i && tree.right[p] != i
            error("Parent-child mismatch: node $i names parent $p, but parent $p does not list node $i as a child.")
        end
    end

    return nothing
end

_root_nodes(tree::Tree) = [i for i in eachindex(tree) if tree.parent[i] == 0]

function _validate_unique_children(tree::Tree)
    child_parent = Dict{Int, Int}()

    for i in eachindex(tree)
        l = tree.left[i]
        r = tree.right[i]

        l != 0 && r != 0 && l == r && error("Node $i cannot list child $l in both left and right.")

        for child in (l, r)
            child == 0 && continue
            if haskey(child_parent, child)
                error("Node $child is listed as a child of both node $(child_parent[child]) and node $i.")
            end
            child_parent[child] = i
        end
    end

    return nothing
end

function _validate_roots(tree::Tree, roots::Vector{Int}; require_single_root::Bool)
    if !isempty(tree) && isempty(roots)
        error("Tree must have at least one root node.")
    end

    if require_single_root && length(roots) != 1
        error("Tree must have exactly one root node; found $(length(roots)).")
    end

    for i in eachindex(tree)
        if tree.parent[i] == 0
            tree.kind[i] == Root || error("Node $i has no parent but kind $(tree.kind[i]); root nodes must have kind Root.")
        elseif tree.kind[i] == Root
            error("Node $i has kind Root but parent $(tree.parent[i]); root nodes are identified by parent == 0.")
        end
    end

    return nothing
end

function _validate_degrees(tree::Tree)
    for i in eachindex(tree)
        k = tree.kind[i]
        p = tree.parent[i]
        nchildren = length(children(tree, i))

        if k == Root
            p == 0 || error("Root node $i must not have a parent.")
            nchildren >= 1 || error("Root node $i must have at least one child.")
        elseif k == Binary
            p != 0 || error("Binary node $i must have a parent.")
            nchildren == 2 || error("Binary node $i must have exactly two children.")
        elseif k == UnsampledUnary || k == SampledUnary
            p != 0 || error("Unary node $i must have a parent.")
            nchildren == 1 || error("Unary node $i must have exactly one child.")
        elseif k == SampledLeaf
            p != 0 || error("SampledLeaf node $i must have a parent.")
            nchildren == 0 || error("SampledLeaf node $i must not have children.")
        else
            error("Unknown NodeKind at node $i.")
        end
    end
    return nothing
end

function _validate_reachable(tree::Tree, roots::Vector{Int})
    seen = Set{Int}()
    stack = reverse(roots)
    while !isempty(stack)
        current = pop!(stack)
        current in seen && continue
        push!(seen, current)
        append!(stack, reverse(children(tree, current)))
    end

    length(seen) == length(tree) || error("Tree contains nodes that are not reachable from root nodes.")
    return nothing
end

function _validate_time_ordering(tree::Tree)
    for i in 2:length(tree)
        tree.time[i - 1] <= tree.time[i] || error("Tree nodes must be stored in nondecreasing time order; node $(i - 1) has time $(tree.time[i - 1]) and node $i has time $(tree.time[i]).")
    end

    for i in eachindex(tree)
        for child in children(tree, i)
            i < child || error("Parent index must be less than child index for canonical time-ordered storage; edge $i -> $child violates this.")
            tree.time[i] < tree.time[child] || error("Parent time must be less than child time for edge $i -> $child.")
        end
    end
    return nothing
end

function _validate_acyclic(tree::Tree)
    for i in eachindex(tree)
        seen = Set{Int}()
        current = i
        while current != 0
            current in seen && error("Cycle detected starting at node $i.")
            push!(seen, current)
            current = tree.parent[current]
        end
    end
    return nothing
end
