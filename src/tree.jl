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

Canonical time-ordered rooted tree storage.

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

function Base.show(io::IO, tree::Tree)
    nroots = count(==(0), tree.parent)
    ninternal = count(i -> tree.left[i] != 0 || tree.right[i] != 0, eachindex(tree))
    ntips = count(i -> tree.left[i] == 0 && tree.right[i] == 0, eachindex(tree))
    print(io,
          "Tree(",
          length(tree),
          " nodes, ",
          nroots,
          " root",
          nroots == 1 ? "" : "s",
          ", ",
          ninternal,
          " internal, ",
          ntips,
          " tips)")
end

function Base.show(io::IO, node::Node)
    print(io, "Node(", node.id, ", ", node.kind, ", time=", node.time)
    node.parent != 0 && print(io, ", parent=", node.parent)

    children = Int[]
    node.left != 0 && push!(children, node.left)
    node.right != 0 && push!(children, node.right)
    !isempty(children) && print(io, ", children=", children)

    node.host != 0 && print(io, ", host=", node.host)
    node.label != 0 && print(io, ", label=", node.label)

    print(io, ")")
end

struct PreorderIterator
    tree::Tree
    starts::Vector{Int}
end

struct PostorderIterator
    tree::Tree
    starts::Vector{Int}
end

struct BreadthFirstIterator
    tree::Tree
    starts::Vector{Int}
end

"""
    Tree()

Construct an empty `Tree`.
"""
Tree() = Tree(Float64[], Int[], Int[], Int[], NodeKind[], Int[], Int[])

Base.length(tree::Tree) = length(tree.time)
Base.isempty(tree::Tree) = isempty(tree.time)
Base.size(tree::Tree) = (length(tree),)
Base.eachindex(tree::Tree) = eachindex(tree.time)
Base.firstindex(::Tree) = 1
Base.lastindex(tree::Tree) = length(tree)
Base.eltype(::Type{Tree}) = Node
Base.IteratorSize(::Type{Tree}) = Base.HasLength()
Base.IteratorEltype(::Type{Tree}) = Base.HasEltype()

"""
    tree[i]

Return a `Node` row view for node id `i`.

The returned `Node` is immutable and mirrors the stored fields for that node.
Node ids are the vector indices in `tree`.
"""
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

"""
    children(tree::Tree, i::Int)

Return the child node ids of node `i` in left-before-right order.

Leaf or tip nodes return `Int[]`.
"""
function children(tree::Tree, i::Int)
    _check_node(tree, i)
    out = Int[]
    tree.left[i] != 0 && push!(out, tree.left[i])
    tree.right[i] != 0 && push!(out, tree.right[i])
    return out
end

"""
    roots(tree::Tree)

Return node ids with no parent.

Canonical validated trees contain exactly one root. Provisional structures may
contain zero or multiple roots.
"""
function roots(tree::Tree)
    return [i for i in eachindex(tree) if tree.parent[i] == 0]
end

"""
    root(tree::Tree)

Return the unique root node id.

Throws if `tree` does not have exactly one root.
"""
function root(tree::Tree)
    rs = roots(tree)
    length(rs) == 1 || error("Tree must have exactly one root; found $(length(rs)).")
    return only(rs)
end

"""
    parent(tree::Tree, i::Int)

Return the parent node id of node `i`, or `0` if `i` is a root.
"""
function parent(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.parent[i]
end

"""
    isleaf(tree::Tree, i::Int)

Return `true` when node `i` has no children.

In TreeSim, a leaf is the structural condition of having no children. Sampled
tips in canonical trees use `SampledLeaf`, but the predicate itself checks the
stored child links.
"""
function isleaf(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.left[i] == 0 && tree.right[i] == 0
end

"""
    isinternal(tree::Tree, i::Int)

Return `true` when node `i` has at least one child.
"""
isinternal(tree::Tree, i::Int) = !isleaf(tree, i)

"""
    isbinary(tree::Tree, i::Int)

Return `true` when node `i` has both a left and right child.
"""
function isbinary(tree::Tree, i::Int)
    _check_node(tree, i)
    return tree.left[i] != 0 && tree.right[i] != 0
end

"""
    isunary(tree::Tree, i::Int)

Return `true` when node `i` has exactly one child.
"""
function isunary(tree::Tree, i::Int)
    _check_node(tree, i)
    return (tree.left[i] != 0) ⊻ (tree.right[i] != 0)
end

"""
    leaves(tree::Tree)

Return node ids with no children, ordered by node id.

For preorder tip encounter order, use `tips(tree)`.
"""
function leaves(tree::Tree)
    return [i for i in eachindex(tree) if isleaf(tree, i)]
end

"""
    tips(tree::Tree)

Return leaf/tip node ids in preorder tip encounter order.

Empty trees return `Int[]`. Multi-root provisional structures are traversed in
`roots(tree)` order.
"""
function tips(tree::Tree)
    out = Int[]
    for i in preorder(tree)
        isleaf(tree, i) && push!(out, i)
    end
    return out
end

"""
    internal_nodes(tree::Tree)

Return node ids with at least one child, ordered by node id.
"""
function internal_nodes(tree::Tree)
    return [i for i in eachindex(tree) if isinternal(tree, i)]
end

"""
    nnodes(tree::Tree)

Return the number of nodes stored in `tree`.
"""
nnodes(tree::Tree) = length(tree)

"""
    nleaves(tree::Tree)

Return the number of leaf/tip nodes, defined structurally as nodes with no
children.
"""
function nleaves(tree::Tree)
    count = 0
    for i in eachindex(tree)
        count += (tree.left[i] == 0 && tree.right[i] == 0)
    end
    return count
end

"""
    ninternal(tree::Tree)

Return the number of internal nodes, defined structurally as nodes with at
least one child. Roots with children count as internal nodes.
"""
function ninternal(tree::Tree)
    count = 0
    for i in eachindex(tree)
        count += (tree.left[i] != 0 || tree.right[i] != 0)
    end
    return count
end

"""
    preorder(tree::Tree[, start])

Iterate over node ids in depth-first preorder from `start`, visiting the left
child before the right child. If `start` is omitted, traversal starts from all
roots in `roots(tree)` order. Empty trees yield no nodes.
"""
preorder(tree::Tree) = PreorderIterator(tree, roots(tree))

function preorder(tree::Tree, start::Int)
    _check_node(tree, start)
    return PreorderIterator(tree, [start])
end

"""
    postorder(tree::Tree[, start])

Iterate over node ids in depth-first postorder from `start`, visiting children
before their parent. Left children are visited before right children. If
`start` is omitted, traversal starts from all roots in `roots(tree)` order.
Empty trees yield no nodes.
"""
postorder(tree::Tree) = PostorderIterator(tree, roots(tree))

function postorder(tree::Tree, start::Int)
    _check_node(tree, start)
    return PostorderIterator(tree, [start])
end

"""
    breadthfirst(tree::Tree[, start])

Iterate over node ids level by level from `start`, visiting the left child
before the right child. If `start` is omitted, traversal starts from all roots
in `roots(tree)` order. Empty trees yield no nodes.
"""
breadthfirst(tree::Tree) = BreadthFirstIterator(tree, roots(tree))

function breadthfirst(tree::Tree, start::Int)
    _check_node(tree, start)
    return BreadthFirstIterator(tree, [start])
end

Base.IteratorSize(::Type{PreorderIterator}) = Base.SizeUnknown()
Base.IteratorSize(::Type{PostorderIterator}) = Base.SizeUnknown()
Base.IteratorSize(::Type{BreadthFirstIterator}) = Base.SizeUnknown()
Base.eltype(::Type{PreorderIterator}) = Int
Base.eltype(::Type{PostorderIterator}) = Int
Base.eltype(::Type{BreadthFirstIterator}) = Int

function Base.iterate(iter::PreorderIterator)
    stack = reverse(copy(iter.starts))
    return iterate(iter, stack)
end

function Base.iterate(iter::PreorderIterator, stack::Vector{Int})
    isempty(stack) && return nothing
    current = pop!(stack)
    right = iter.tree.right[current]
    left = iter.tree.left[current]
    right != 0 && push!(stack, right)
    left != 0 && push!(stack, left)
    return (current, stack)
end

function Base.iterate(iter::PostorderIterator)
    stack = Tuple{Int, Bool}[]
    for start in reverse(iter.starts)
        push!(stack, (start, false))
    end
    return iterate(iter, stack)
end

function Base.iterate(iter::PostorderIterator, stack::Vector{Tuple{Int, Bool}})
    while !isempty(stack)
        current, expanded = pop!(stack)
        if expanded
            return (current, stack)
        end

        push!(stack, (current, true))
        right = iter.tree.right[current]
        left = iter.tree.left[current]
        right != 0 && push!(stack, (right, false))
        left != 0 && push!(stack, (left, false))
    end
    return nothing
end

function Base.iterate(iter::BreadthFirstIterator)
    return iterate(iter, (copy(iter.starts), 1))
end

function Base.iterate(iter::BreadthFirstIterator, state::Tuple{Vector{Int}, Int})
    queue, offset = state
    offset > length(queue) && return nothing
    current = queue[offset]
    left = iter.tree.left[current]
    right = iter.tree.right[current]
    left != 0 && push!(queue, left)
    right != 0 && push!(queue, right)
    return (current, (queue, offset + 1))
end

"""
    ancestors(tree::Tree, i::Int)

Return parent ids on the path from node `i` toward the root.

The immediate parent comes first. Roots return `Int[]`.
"""
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

"""
    descendants(tree::Tree, i::Int)

Return descendants of node `i` in preorder, excluding `i` itself.
"""
function descendants(tree::Tree, i::Int)
    _check_node(tree, i)
    out = collect(preorder(tree, i))
    return out[2:end]
end

"""
    subtree_nodes(tree::Tree, i::Int)

Return all node ids in the subtree rooted at `i`, ordered by preorder and
including `i` itself.
"""
function subtree_nodes(tree::Tree, i::Int)
    return collect(preorder(tree, i))
end

"""
    node_depths(tree::Tree)

Return a node-index-aligned vector of metric depths. `depths[i]` is the stored
branch-length distance from node `i`'s own root to node `i`. Empty trees return
`Float64[]`; multi-root structures are measured within each rooted component.
"""
function node_depths(tree::Tree)
    depths = zeros(Float64, length(tree))

    for i in preorder(tree)
        p = tree.parent[i]
        depths[i] = p == 0 ? 0.0 : depths[p] + tree.time[i] - tree.time[p]
    end

    return depths
end

"""
    root_to_tip_distances(tree::Tree)

Return root-to-tip distances for tips only.

Distances are ordered by preorder tip encounter order, matching `tips(tree)`.
In multi-root provisional structures, each distance is measured from the tip's
own root. Empty trees return `Float64[]`.
"""
function root_to_tip_distances(tree::Tree)
    depths = zeros(Float64, length(tree))
    out = Float64[]

    for i in preorder(tree)
        p = tree.parent[i]
        depths[i] = p == 0 ? 0.0 : depths[p] + tree.time[i] - tree.time[p]
        if tree.left[i] == 0 && tree.right[i] == 0
            push!(out, depths[i])
        end
    end

    return out
end

"""
    tree_height(tree::Tree)

Return the maximum root-to-tip distance. Empty trees have height `0.0`. For
multi-root structures, the height is the maximum across rooted components.
"""
function tree_height(tree::Tree)
    height = 0.0
    for distance in root_to_tip_distances(tree)
        height = max(height, distance)
    end
    return height
end

"""
    mean_root_to_tip_distance(tree::Tree)

Return the mean root-to-tip distance across all tips in preorder tip encounter
order. Empty trees return `0.0`; multi-root structures include tips from all
rooted components.
"""
function mean_root_to_tip_distance(tree::Tree)
    distances = root_to_tip_distances(tree)
    isempty(distances) && return 0.0
    return sum(distances) / length(distances)
end

"""
    ncherries(tree::Tree)

Return the number of cherries. A cherry is an internal node whose left and right
children are both tips; unary nodes do not count as cherries.
"""
function ncherries(tree::Tree)
    count = 0
    for i in eachindex(tree)
        l = tree.left[i]
        r = tree.right[i]
        count += (l != 0 && r != 0 && isleaf(tree, l) && isleaf(tree, r))
    end
    return count
end

"""
    tip_positions(tree::Tree)

Return a node-index-aligned vector of vertical tip positions. Tip entries are
consecutive `Float64` values in preorder tip encounter order, starting at
`1.0`; non-tip entries are `NaN`. Empty trees return `Float64[]`.
"""
function tip_positions(tree::Tree)
    y = fill(NaN, length(tree))

    position = 1.0
    for i in tips(tree)
        y[i] = position
        position += 1.0
    end

    return y
end

"""
    node_positions(tree::Tree)

Return `(x, y)` node-index-aligned coordinate vectors for derived layouts.
`x == node_depths(tree)`. Tip `y` values come from `tip_positions(tree)`;
internal node `y` values are the mean of their child `y` values. Empty trees
return `(Float64[], Float64[])`.
"""
function node_positions(tree::Tree)
    x = node_depths(tree)
    y = tip_positions(tree)

    for i in postorder(tree)
        isleaf(tree, i) && continue

        total = 0.0
        nchildren = 0
        l = tree.left[i]
        r = tree.right[i]

        if l != 0
            total += y[l]
            nchildren += 1
        end
        if r != 0
            total += y[r]
            nchildren += 1
        end

        nchildren > 0 && (y[i] = total / nchildren)
    end

    return x, y
end

"""
    parent_child_pairs(tree::Tree)

Return directed tree edges as `(parent, child)` vectors. Edges are ordered by
parent node index, with each parent's left child emitted before its right child.
Absent child slots are skipped. Empty trees return `(Int[], Int[])`.
"""
function parent_child_pairs(tree::Tree)
    parents = Int[]
    children = Int[]

    for i in eachindex(tree)
        l = tree.left[i]
        r = tree.right[i]

        if l != 0
            push!(parents, i)
            push!(children, l)
        end
        if r != 0
            push!(parents, i)
            push!(children, r)
        end
    end

    return parents, children
end

"""
    edge_segments(tree::Tree)

Return `(x1, y1, x2, y2)` coordinate vectors for each directed edge. Segment
order exactly matches `parent_child_pairs(tree)`, and coordinates are derived
from `node_positions(tree)`. Empty trees return four empty `Float64` vectors.
"""
function edge_segments(tree::Tree)
    parents, children = parent_child_pairs(tree)
    x, y = node_positions(tree)

    x1 = Float64[]
    y1 = Float64[]
    x2 = Float64[]
    y2 = Float64[]

    for k in eachindex(parents)
        p = parents[k]
        c = children[k]
        push!(x1, x[p])
        push!(y1, y[p])
        push!(x2, x[c])
        push!(y2, y[c])
    end

    return x1, y1, x2, y2
end

"""
    branch_length(tree::Tree, i::Int)

Return the stored branch length from `parent(tree, i)` to node `i`.

Throws for a root node because roots have no parent branch.
"""
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
