module TreeSim

include("tree.jl")
include("plot.jl")

"""
    tree_from_eventlog(log; validate=true)

Convert a supported epidemic event log into one sampled-ancestry `Tree`.

This package defines the API; concrete event-log adaptors are loaded by package
extensions so core TreeSim remains independent of epidemic simulation packages.

This is the strict single-tree convenience API. It returns `Tree()` when no
sampled ancestry is retained, returns the only retained component when exactly
one exists, and errors if the log retains multiple seeded sampled components.
Use `forest_from_eventlog` to extract all retained sampled components
explicitly.

The bridge extracts sampled transmission ancestry, not complete outbreak
history: unobserved side lineages without sampled descendants are pruned.
Retained parent-child event times must be strictly increasing because canonical
TreeSim edges have positive temporal length; equal-time retained ancestry is
reported as an error rather than silently adjusted.
"""
function tree_from_eventlog end

"""
    forest_from_eventlog(log; validate=true)

Convert a supported epidemic event log into all retained sampled ancestry
components as a `Vector{Tree}`.

This is the full sampled-ancestry extraction API. It returns `Tree[]` when no
sampled ancestry is retained, otherwise one standalone canonical `Tree` per
retained seeded sampled component. For EpiSim event logs, components are sorted
deterministically by root time, then root host id, then original root index.

As with `tree_from_eventlog`, this is not a full outbreak-history
reconstruction: unobserved side lineages without sampled descendants are
pruned, and retained ancestry must satisfy TreeSim's strict temporal ordering.
"""
function forest_from_eventlog end

"""
    validate_tree_against_eventlog(log, tree_or_forest)

Validate that a sampled-ancestry `Tree` or `Vector{Tree}` preserves the
bridge-specific correspondence to a supported epidemic event log.

For forests, validation uses one global event-count ledger across all
components so that an event cannot be consumed by more than one component.
Validation also checks TreeSim structure, sampled-event correspondence,
transmission-event correspondence, pruned-node expectations, and strict
retained ancestry timing.
"""
function validate_tree_against_eventlog end

export # Core tree types and node roles
       NodeKind,
       Root,
       Binary,
       UnsampledUnary,
       SampledUnary,
       SampledLeaf,
       Node,
       Tree,

       # Tree and node inspection
       root,
       roots,
       children,
       parent,
       isleaf,
       isinternal,
       isbinary,
       isunary,
       leaves,
       tips,
       internal_nodes,
       nnodes,
       nleaves,
       ninternal,

       # Traversal
       preorder,
       postorder,
       breadthfirst,
       ancestors,
       descendants,
       subtree_nodes,

       # Basic summaries
       branch_length,
       node_depths,
       root_to_tip_distances,
       tree_height,
       mean_root_to_tip_distance,
       ncherries,

       # Layout and plotting
       tip_positions,
       node_positions,
       parent_child_pairs,
       edge_segments,
       plot_tree,

       # Validation
       validate_tree,

       # Optional event-log adaptors
       tree_from_eventlog,
       forest_from_eventlog,
       validate_tree_against_eventlog

end
