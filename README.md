# TreeSim.jl

`TreeSim.jl` is the canonical standalone tree representation package for the
recovered outbreak-modelling ecosystem.

The current package scope is deliberately small. It provides:

- a canonical `Tree` type
- structural node roles via `NodeKind`
- basic root, parent, child, ancestor, descendant, tip, and leaf utilities
- preorder, postorder, and breadth-first traversal iterators
- lightweight node-count and root-to-tip distance statistics
- validation for the structural invariants downstream packages rely on

It does not provide birth-death likelihoods, epidemic simulation, sequence
simulation, plotting, large import/export systems, or orchestration logic.

## Canonical Tree Contract

A `Tree` stores one node per vector index. Node ids are the vector indices, and
`0` is the null sentinel for absent parent or child references.

The active canonical contract is strict:

- non-empty validated trees have exactly one `Root`
- every node is reachable from the root
- parent and child references are reciprocal
- a child may appear in only one parent slot
- node degree must match `NodeKind`
- nodes are stored in nondecreasing time order
- every edge points from a lower index to a higher index
- parent time is strictly less than child time

Time-ordered storage is a hard invariant of the canonical representation, not a
formatting preference.

## Node Kinds

- `Root`: unique node with `parent == 0`; it must have at least one child.
- `Binary`: non-root internal node with exactly two children.
- `UnsampledUnary`: non-root unsampled internal node with exactly one child.
- `SampledUnary`: non-root sampled internal node with exactly one child.
- `SampledLeaf`: non-root sampled tip with no children.

`SampledUnary` and `SampledLeaf` encode structural sampling status only. This
package does not decide whether a tree is admissible for any particular
birth-death, epidemiological, or sequence model.

## Example

```julia
using TreeSim

tree = Tree(
    [0.0, 0.4, 0.8],
    [2, 0, 0],
    [3, 0, 0],
    [0, 1, 1],
    [Root, SampledLeaf, SampledLeaf],
    [0, 0, 0],
    [0, 101, 102],
)

validate_tree(tree)
children(tree, 1)      # [2, 3]
branch_length(tree, 2) # 0.4
```

## Traversal Helpers

Traversal functions yield node ids and preserve the canonical left-before-right
child order. `root(tree)` is strict and returns only a unique root; it throws
for empty or multi-root trees. `roots(tree)` returns all roots in node-index
order.

```julia
collect(preorder(tree))     # [1, 2, 3]
collect(postorder(tree))    # [2, 3, 1]
collect(breadthfirst(tree)) # [1, 2, 3]

root(tree)           # 1
roots(tree)          # [1]
leaves(tree)         # [2, 3]
tips(tree)           # [2, 3]
internal_nodes(tree) # [1]
```

`preorder(tree)`, `postorder(tree)`, and `breadthfirst(tree)` traverse all roots
in `roots(tree)` order when no start node is given. Empty trees yield empty
traversals. Pass an explicit start node, such as `preorder(tree, 2)`, to
traverse a subtree or to inspect provisional structures before full validation.
`descendants(tree, i)` is an eager convenience collector using the same
left-before-right preorder contract, excluding `i` itself.

`tips(tree)` is an eager convenience collector for leaf nodes in preorder tip
encounter order. This is the tip-ordering convention used by downstream
analysis helpers. `subtree_nodes(tree, i)` eagerly collects the subtree rooted
at `i` in preorder and includes `i` itself; it does not rebuild or copy a tree.

## Basic Statistics

The first statistics layer is intentionally small and derives from the core
tree representation:

```julia
nnodes(tree)               # length(tree)
nleaves(tree)              # tips with no children
ninternal(tree)            # nodes with at least one child, including roots
node_depths(tree)          # metric depth by node index
root_to_tip_distances(tree) # tips in preorder encounter order
tree_height(tree)          # maximum root-to-tip distance
mean_root_to_tip_distance(tree)
ncherries(tree)
```

`root_to_tip_distances(tree)` returns distances for tips only. Tips are ordered
by preorder encounter, and in multi-root structures each distance is measured
from that tip's own root. `tree_height(tree)` is the maximum of those distances;
empty trees have height `0.0`.

`node_depths(tree)` returns a node-index-aligned vector of metric depths using
stored branch lengths, so `node_depths(tree)[i]` is the distance from node `i`'s
own root to node `i`. Empty trees return `Float64[]`, and multi-root structures
are measured within each rooted component.

`mean_root_to_tip_distance(tree)` is the mean over all tip distances and returns
`0.0` for empty trees. `ncherries(tree)` counts internal nodes whose left and
right children are both tips; unary nodes do not count as cherries.

## Validation Notes

`validate_tree(tree)` validates the canonical trusted representation. The
keyword arguments `require_single_root=false` and `require_reachable=false` are
available for inspecting provisional non-canonical structures, but downstream
packages should rely on the default strict validation path.

See `VALIDATION.md` and `TRUST_CRITERIA.md` for the package recovery criteria
guiding this minimal core.
