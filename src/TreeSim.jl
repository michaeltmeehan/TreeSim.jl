module TreeSim

include("tree.jl")
include("plot.jl")

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
       validate_tree

end
