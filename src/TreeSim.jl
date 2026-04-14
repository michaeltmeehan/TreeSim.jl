module TreeSim

include("tree.jl")

export NodeKind,
       Root,
       Binary,
       UnsampledUnary,
       SampledUnary,
       SampledLeaf,
       Node,
       Tree,
       children,
       parent,
       isleaf,
       isinternal,
       isbinary,
       isunary,
       leaves,
       ancestors,
       descendants,
       branch_length,
       validate_tree

end
