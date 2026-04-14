# TreeSim.jl

`TreeSim.jl` is the canonical standalone tree representation package for the
recovered outbreak-modelling ecosystem.

The current package scope is deliberately small. It provides:

- a canonical `Tree` type
- structural node roles via `NodeKind`
- basic parent, child, ancestor, descendant, and leaf utilities
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

## Validation Notes

`validate_tree(tree)` validates the canonical trusted representation. The
keyword arguments `require_single_root=false` and `require_reachable=false` are
available for inspecting provisional non-canonical structures, but downstream
packages should rely on the default strict validation path.

See `VALIDATION.md` and `TRUST_CRITERIA.md` for the package recovery criteria
guiding this minimal core.
