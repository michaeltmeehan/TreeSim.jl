module TreeSim

include("tree.jl")
include("statistics.jl")
include("plot.jl")

"""
    tree_from_eventlog(log; validate=true)

Convert a supported epidemic event log into one reconstructed `Tree`.

This package defines the API; concrete event-log adaptors are loaded by package
extensions so core TreeSim remains independent of epidemic simulation packages.

This is the strict single-tree convenience API. It returns `Tree()` when no
reconstructed ancestry is retained, returns the only retained component when exactly
one exists, and errors if the log retains multiple seeded sampled components.
Use `forest_from_eventlog` to extract all retained sampled components
explicitly.

`tree_from_eventlog` is an alias for `reconstructed_tree_from_eventlog`.
Reconstructed extraction keeps sampled tips and sampled internal ancestry, but
collapses unsampled unary transmission intermediates that belong to the
generative outbreak history rather than the reconstructed tree used for
inference.
"""
function tree_from_eventlog end

"""
    reconstructed_tree_from_eventlog(log; validate=true)

Convert a supported epidemic event log into one reconstructed sampled tree.

Reconstructed trees retain serial/fossilized sampling observations and the
branching/sampled internal ancestry needed to connect them. Unsampled unary
transmission intermediates are collapsed at extraction time.
"""
function reconstructed_tree_from_eventlog end

"""
    full_tree_from_eventlog(log; validate=true)

Convert a supported epidemic event log into one fuller sampled-ancestry tree.

This extractor preserves unsampled unary transmission intermediates that are
ancestral to retained sampled nodes. It is still represented with TreeSim's
current sampled-tree node kinds, so completely unsampled removal/activation
side histories are not materialized as tips.
"""
function full_tree_from_eventlog end

"""
    forest_from_eventlog(log; validate=true)

Convert a supported epidemic event log into all reconstructed ancestry
components as a `Vector{Tree}`.

This is the reconstructed forest extraction API. It returns `Tree[]` when no
sampled ancestry is retained, otherwise one standalone canonical `Tree` per
retained seeded sampled component. For EpiSim event logs, components are sorted
deterministically by root time, then root host id, then original root index.

As with `tree_from_eventlog`, this is not a full outbreak-history
reconstruction.
"""
function forest_from_eventlog end

"""
    validate_tree_against_eventlog(log, tree_or_forest)

Validate that a reconstructed `Tree` or `Vector{Tree}` preserves the
bridge-specific correspondence to a supported epidemic event log.

For forests, validation uses one global event-count ledger across all
components so that an event cannot be consumed by more than one component.
Validation also checks TreeSim structure, sampled-event correspondence,
reconstructed transmission-event correspondence, pruned-node expectations, and
strict retained ancestry timing.
"""
function validate_tree_against_eventlog end

"""
    validate_full_tree_against_eventlog(log, tree_or_forest)

Validate that a fuller sampled-ancestry `Tree` or `Vector{Tree}` preserves the
bridge-specific correspondence to a supported epidemic event log.

Full-tree validation preserves the stricter transmission correspondence used by
`full_tree_from_eventlog`: retained binary and unsampled-unary transmission
nodes must match immediate event-log infectee/infector/time signatures.
"""
function validate_full_tree_against_eventlog end

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
       TreeStats,
       tree_statistics,

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
       reconstructed_tree_from_eventlog,
       full_tree_from_eventlog,
       forest_from_eventlog,
       validate_tree_against_eventlog,
       validate_full_tree_against_eventlog

end
