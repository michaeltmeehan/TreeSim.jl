# TreeSim.jl

`TreeSim.jl` is a small Julia package for representing, validating, traversing,
summarizing, and quickly viewing rooted trees.

It is a canonical tree layer for outbreak-modelling workflows. Its job is
deliberately narrow: keep tree storage efficient and unambiguous, provide
reliable traversal and structural checks, and expose modest analytical helpers
that downstream packages can build on.

`TreeSim.jl` does not simulate birth-death trees, run epidemic models, simulate
sequences, or estimate likelihoods.

## Installation

From the package repository:

```julia
using Pkg
Pkg.activate("repos/TreeSim.jl")
Pkg.instantiate()
Pkg.test()
```

From another local project that should use this checkout:

```julia
using Pkg
Pkg.develop(path="repos/TreeSim.jl")
using TreeSim
```

Adjust the path if your Julia session starts inside `repos/TreeSim.jl`; in that
case use `Pkg.activate(".")` or `Pkg.develop(path=".")`.

## Core Concepts

A `Tree` stores one node per vector index. Node ids are the vector indices, and
`0` is the null value for an absent parent or child.

Each node has:

- `time`: node time
- `left` and `right`: child node ids, or `0`
- `parent`: parent node id, or `0` for a root
- `kind`: a `NodeKind` describing the structural role
- `host` and `label`: integer metadata fields for downstream packages

The available node kinds are:

- `Root`: a root node with at least one child
- `Binary`: a non-root internal node with exactly two children
- `UnsampledUnary`: an unsampled non-root internal node with one child
- `SampledUnary`: a sampled non-root internal node with one child
- `SampledLeaf`: a sampled non-root tip with no children

Validated canonical trees are strict: they have one root, all nodes are
reachable from that root, parent and child references agree, node degrees match
`NodeKind`, node times are stored in nondecreasing order, and every edge points
from a lower node index to a higher node index with strictly increasing time.

## Worked Example

The examples below all use the same five-node tree. It has one root, one
internal branching node, and three sampled tips.

```julia
using TreeSim

tree = Tree(
    [0.0, 0.4, 0.8, 1.1, 1.3],                 # time
    [2,   4,   0,   0,   0],                   # left child
    [3,   5,   0,   0,   0],                   # right child
    [0,   1,   1,   2,   2],                   # parent
    [Root, Binary, SampledLeaf, SampledLeaf, SampledLeaf],
    [0, 0, 0, 0, 0],                           # host metadata
    [0, 0, 101, 102, 103],                     # label metadata
)
```

Inspecting `tree` gives a compact structural summary:

```julia
tree
```

```text
Tree(5 nodes, 1 root, 2 internal, 3 tips)
```

The summary is intentionally small: it tells you the stored size, root count,
internal-node count, and tip count without hiding the efficient vector-backed
representation.

## Validate The Structure

Use `validate_tree` before treating a hand-built or imported tree as canonical.

```julia
validate_tree(tree) # true
```

Validation checks the structural contract that downstream packages rely on:
there is one root, every node is reachable, parent and child references agree,
node kinds match their degree, storage is time ordered, and edge times strictly
increase away from the root.

## Inspect Nodes And Relationships

Indexing a tree returns an immutable `Node` row view. The display emphasizes the
meaningful relationship fields for that node kind.

```julia
tree[1] # Node(1, Root, time=0.0, children=[2, 3])
tree[2] # Node(2, Binary, time=0.4, parent=1, children=[4, 5])
tree[4] # Node(4, SampledLeaf, time=1.1, parent=2)
```

Relationship helpers return node ids:

```julia
root(tree)            # 1
children(tree, 2)     # [4, 5]
parent(tree, 4)       # 2
ancestors(tree, 5)    # [2, 1]
descendants(tree, 2)  # [4, 5]
subtree_nodes(tree, 2) # [2, 4, 5]
```

Node `2` is the internal branching node. Its descendants are tips `4` and `5`,
and `subtree_nodes(tree, 2)` includes node `2` itself because it names the
entire subtree rooted there.

## Traverse The Tree

Traversal functions yield node ids and visit left children before right
children.

```julia
collect(preorder(tree))      # [1, 2, 4, 5, 3]
collect(postorder(tree))     # [4, 5, 2, 3, 1]
collect(breadthfirst(tree))  # [1, 2, 3, 4, 5]
```

Preorder is useful when you want to encounter a parent before its descendants.
Postorder is useful for bottom-up summaries. Breadth-first traversal visits
nodes level by level.

You can also start traversal from a subtree:

```julia
collect(preorder(tree, 2))  # [2, 4, 5]
```

Tip helpers use preorder tip encounter order:

```julia
tips(tree)           # [4, 5, 3]
internal_nodes(tree) # [1, 2]
```

This tip ordering is also used by root-to-tip statistics and layout helpers.

## Summarize The Tree

The statistics helpers are deliberately modest and derive directly from the
stored tree.

```julia
nnodes(tree)    # 5
ninternal(tree) # 2
nleaves(tree)   # 3
ncherries(tree) # 1
```

`ncherries(tree)` returns `1` because node `2` has two tip children. The root
does not count as a cherry because one of its children, node `2`, is internal.

Depth and branch-length helpers use stored node times:

```julia
branch_length(tree, 4)          # 0.7000000000000001
node_depths(tree)               # [0.0, 0.4, 0.8, 1.1, 1.3000000000000003]
root_to_tip_distances(tree)     # [1.1, 1.3000000000000003, 0.8]
tree_height(tree)               # 1.3000000000000003
mean_root_to_tip_distance(tree) # 1.0666666666666667
```

`node_depths(tree)[4]` is `1.1`, the distance from the root to node `4`.
`root_to_tip_distances(tree)` reports tips in the same order as `tips(tree)`,
so the distances correspond to tips `[4, 5, 3]`.

## Plot For Quick Inspection

`plot_tree(tree)` returns a lightweight SVG-backed `TreePlot` display object for
quick inspection in notebooks, VS Code, Pluto, or any display that can render
SVG.

```julia
plot_tree(tree)
plot_tree(tree; labels=string.(eachindex(tree)))
plot_tree(tree; node_color=tree.host)
plot_tree(tree; width=640, height=360, edge_color="black")
```

The viewer uses the same traversal and layout conventions described above.
`node_color` and `labels` must be node-index-aligned when supplied. Numeric
colors use a continuous scale; other values use discrete categories.

This is not a plotting framework. Per-edge aesthetics, interactive plotting,
layout redesign, file export helpers, and plotting-framework integrations are
outside the current package scope.

## Reuse Layout Coordinates

When downstream code needs coordinates without a display object, use the layout
helpers directly.

```julia
tip_positions(tree)       # [NaN, NaN, 3.0, 1.0, 2.0]
node_positions(tree)      # (x, y), with x == node_depths(tree)
parent_child_pairs(tree)  # ([1, 1, 2, 2], [2, 3, 4, 5])
edge_segments(tree)       # (x1, y1, x2, y2)
```

`tip_positions(tree)` assigns y-values to tips in preorder tip order, so tips
`[4, 5, 3]` receive positions `[1.0, 2.0, 3.0]`. Internal node positions are the
mean of their child positions.

## Extract Sampled Ancestry From EpiSim

TreeSim keeps epidemic simulators optional. When `EpiSim` is loaded alongside
`TreeSim`, a package extension provides a narrow sampled-ancestry bridge from
`EpiSim.EventLog` to TreeSim trees.

Use `tree_from_eventlog` when the retained sampled ancestry is known to have at
most one seeded component. It returns `Tree()` when there are no retained
samples, returns the only retained component when there is one, and errors if
multiple retained sampled components are present.

```julia
using TreeSim
using EpiSim

log = EventLog(
    [0.0, 1.0, 2.0],
    [1, 2, 2],
    [0, 1, 0],
    [EK_Seeding, EK_Transmission, EK_SerialSampling],
)

tree = tree_from_eventlog(log)
validate_tree_against_eventlog(log, tree)
```

Use `forest_from_eventlog` when an event log may contain multiple independently
seeded sampled components. It returns one standalone canonical `Tree` per
retained component, ordered by root time, then root host id, then original root
index.

```julia
log = EventLog(
    [0.0, 0.0, 1.0, 1.5],
    [1, 2, 1, 2],
    [0, 0, 0, 0],
    [EK_Seeding, EK_Seeding, EK_SerialSampling, EK_SerialSampling],
)

forest = forest_from_eventlog(log)
validate_tree_against_eventlog(log, forest)
```

This bridge extracts sampled transmission ancestry, not the complete outbreak
history. Unsampled side lineages without sampled descendants are pruned. If
retained ancestry would require an equal-time or backward-time edge, conversion
fails rather than adding hidden time perturbations.

## Main Features

Use this map after the worked example when you know the task but not the helper
name.

| Task | Start with |
|---|---|
| Build or inspect trees | `Tree`, `Node`, `tree[i]`, `children`, `parent` |
| Check tree validity | `validate_tree`, `root`, `roots` |
| Find structural groups | `tips`, `leaves`, `internal_nodes`, `isleaf`, `isinternal`, `isbinary`, `isunary` |
| Traverse node ids | `preorder`, `postorder`, `breadthfirst`, `ancestors`, `descendants`, `subtree_nodes` |
| Compute compact summaries | `nnodes`, `nleaves`, `ninternal`, `branch_length`, `node_depths`, `root_to_tip_distances`, `tree_height`, `mean_root_to_tip_distance`, `ncherries` |
| Make a quick visual check | `plot_tree` |
| Reuse layout coordinates | `tip_positions`, `node_positions`, `parent_child_pairs`, `edge_segments` |
| Extract sampled ancestry from EpiSim logs | `tree_from_eventlog`, `forest_from_eventlog`, `validate_tree_against_eventlog` |

`TreeSim.jl` keeps these helpers close to the canonical tree representation.
Simulation, likelihood, and rich plotting workflows belong in downstream
packages built on this tree layer.

## Further Notes

`validate_tree(tree)` defaults to the canonical trusted representation. The
keywords `require_single_root=false` and `require_reachable=false` are available
for inspecting provisional non-canonical structures, but downstream packages
should rely on the default strict validation path.

See `VALIDATION.md` and `TRUST_CRITERIA.md` for the validation contract behind
this minimal core.
