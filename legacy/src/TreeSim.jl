"""
    TreeSim

TreeSim is a Julia package for simulating phylogenetic trees from transmission linelist data. The package provides tools to generate within-host trees, combine them into a single tree, and calculate probabilities and constraints associated with phylogenetic processes.

# Features
- Define and manipulate host structures using the `Host` type.
- Simulate within-host phylogenetic trees.
- Combine multiple within-host trees into a single phylogenetic tree.
- Calculate probabilities and enforce constraints for coalescent events.
- Normalize and relabel phylogenetic trees for analysis.

# Main Functions
- `simulate_phylogeny(linelist::DataFrame; Nₑ::Float64=1e-6)`: Simulate a phylogenetic tree from a linelist DataFrame.
- `get_hosts(linelist::DataFrame)`: Construct a dictionary of `Host` objects from a linelist DataFrame.
"""
module TreeSim

using DataFrames
using DataFramesMeta
using EpiSim
using Lazy
using Phylo
using RecipesBase
using RecipesPipeline
using UnPack


"""
    struct Host

A structure representing a host in a transmission / phylogenetic tree, containing information about the host, the infector, and the timing of various events.

# Fields
- `id::Int64`: The unique identifier for the host.
- `type::Int64`: The type of the host.
- `infector::Union{Int64, Nothing}`: The identifier of the host's infector, or `nothing` if the host has no infector.
- `infector_type::Union{Int64, Nothing}`: The type of the host's infector, or `nothing` if the host has no infector.
- `root::Union{Float64, Nothing}`: The time of the root event, or `nothing` if not applicable.
- `leaf_times::Vector{Float64}`: A vector of times at which the host's leaves were sampled.
- `leaf_ids::Vector{Int64}`: A vector of identifiers for the leaves of the host.
"""
mutable struct Host
	id::Int64
    type::Int64
    infector::Union{Int64, Nothing}
    infector_type::Union{Int64, Nothing}
	root::Union{Float64, Nothing}
	leaf_times::Vector{Float64}
	leaf_ids::Vector{Int64}
end


include("utils.jl")
include("tree_operations.jl")
include("probabilities.jl")
include("constraints.jl")
include("backward.jl")
include("forward.jl")
include("sampling.jl")
include("tree_construction.jl")
include("plot.jl")

export Host, Phylogeny
export simulate_phylogeny, get_hosts


end # module TreeSim
