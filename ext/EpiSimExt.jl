module EpiSimExt

using EpiSim
using TreeSim

const EventLog = EpiSim.EventLog

struct ExtractionSpec
    collapse_unsampled_unary::Bool
end

struct ValidationSpec
    allow_unsampled_unary::Bool
    require_immediate_transmission_child::Bool
end

const RECONSTRUCTED_TREE_SPEC = ExtractionSpec(true)
const FULL_TREE_SPEC = ExtractionSpec(false)
const RECONSTRUCTED_VALIDATION_SPEC = ValidationSpec(false, false)
const FULL_VALIDATION_SPEC = ValidationSpec(true, true)

function TreeSim.tree_from_eventlog(log::EventLog; validate::Bool=true)
    return TreeSim.reconstructed_tree_from_eventlog(log; validate)
end

function TreeSim.reconstructed_tree_from_eventlog(log::EventLog; validate::Bool=true)
    forest = TreeSim.forest_from_eventlog(log; validate)
    isempty(forest) && return TreeSim.Tree()
    length(forest) == 1 ||
        error("reconstructed_tree_from_eventlog found $(length(forest)) retained sampled components. It is a strict single-tree API; use forest_from_eventlog to extract the reconstructed forest.")
    return only(forest)
end

function TreeSim.full_tree_from_eventlog(log::EventLog; validate::Bool=true)
    forest = _forest_from_eventlog(log, FULL_TREE_SPEC; validate)
    isempty(forest) && return TreeSim.Tree()
    length(forest) == 1 ||
        error("full_tree_from_eventlog found $(length(forest)) retained sampled components. It is a strict single-tree API.")
    return only(forest)
end

function TreeSim.forest_from_eventlog(log::EventLog; validate::Bool=true)
    return _forest_from_eventlog(log, RECONSTRUCTED_TREE_SPEC; validate)
end

function _forest_from_eventlog(log::EventLog, spec::ExtractionSpec; validate::Bool=true)
    EpiSim.validate_event_log(log)
    _validate_bridge_log_contract(log)
    forest = _extract_sampled_forest(log, spec)
    validation_spec = spec.collapse_unsampled_unary ? RECONSTRUCTED_VALIDATION_SPEC : FULL_VALIDATION_SPEC
    validate && _validate_forest_against_eventlog(log, forest, validation_spec)
    return forest
end

function TreeSim.validate_tree_against_eventlog(log::EventLog, tree::TreeSim.Tree)
    _validate_single_tree_contract(tree)
    return _validate_forest_against_eventlog(log, isempty(tree) ? TreeSim.Tree[] : [tree], RECONSTRUCTED_VALIDATION_SPEC)
end

function TreeSim.validate_tree_against_eventlog(log::EventLog, forest::AbstractVector{<:TreeSim.Tree})
    return _validate_forest_against_eventlog(log, forest, RECONSTRUCTED_VALIDATION_SPEC)
end

function TreeSim.validate_full_tree_against_eventlog(log::EventLog, tree::TreeSim.Tree)
    _validate_single_tree_contract(tree)
    return _validate_forest_against_eventlog(log, isempty(tree) ? TreeSim.Tree[] : [tree], FULL_VALIDATION_SPEC)
end

function TreeSim.validate_full_tree_against_eventlog(log::EventLog, forest::AbstractVector{<:TreeSim.Tree})
    return _validate_forest_against_eventlog(log, forest, FULL_VALIDATION_SPEC)
end

function _validate_forest_against_eventlog(log::EventLog, forest::AbstractVector{<:TreeSim.Tree}, spec::ValidationSpec)
    EpiSim.validate_event_log(log)
    _validate_bridge_log_contract(log)
    for tree in forest
        _validate_component_contract(tree)
        TreeSim.validate_tree(tree)
    end

    event_counts = Dict{Tuple{EpiSim.EventKind,Int,Float64},Int}()
    transmission_counts = Dict{Tuple{Int,Int,Float64},Int}()
    transmission_by_infector_counts = Dict{Tuple{Int,Float64},Int}()
    sampling_count = 0

    for i in 1:length(log)
        kind = log.kind[i]
        host = log.host[i]
        time = log.time[i]
        _increment!(event_counts, (kind, host, time))

        if kind == EpiSim.EK_Transmission
            _increment!(transmission_counts, (host, log.infector[i], time))
            _increment!(transmission_by_infector_counts, (log.infector[i], time))
        elseif _is_sampling(kind)
            sampling_count += 1
        end
    end

    sampling_count == 0 && isempty(forest) && return true
    sampling_count == 0 && error("Sampled-ancestry tree/forest was provided, but the EventLog has no sampling events.")
    isempty(forest) && error("Sampled-ancestry forest is empty, but the EventLog contains sampling events.")
    seen = Set{Tuple{TreeSim.NodeKind,Int,Float64}}()

    for tree in forest
        for i in eachindex(tree)
            kind = tree.kind[i]
            host = tree.host[i]
            time = tree.time[i]

            if kind == TreeSim.Root
                _consume!(event_counts, (EpiSim.EK_Seeding, host, time)) ||
                    error("Root does not correspond to Seeding event (host=$host, time=$time).")
            elseif kind == TreeSim.SampledLeaf
                (_consume!(event_counts, (EpiSim.EK_SerialSampling, host, time)) ||
                 _consume!(event_counts, (EpiSim.EK_FossilizedSampling, host, time))) ||
                    error("SampledLeaf does not correspond to sampling event (host=$host, time=$time).")
            elseif kind == TreeSim.SampledUnary
                _consume!(event_counts, (EpiSim.EK_FossilizedSampling, host, time)) ||
                    error("SampledUnary does not correspond to FossilizedSampling event (host=$host, time=$time).")
            elseif kind == TreeSim.Binary
                if spec.require_immediate_transmission_child
                    infectees = [tree.host[child] for child in TreeSim.children(tree, i)]
                    _consume_any!(transmission_counts, ((infectee, host, time) for infectee in infectees)) ||
                        error("Binary node at index $i does not match immediate Transmission infectee/infector event.")
                else
                    _consume!(transmission_by_infector_counts, (host, time)) ||
                        error("Binary node at index $i does not match retained Transmission infector/time event.")
                end
            elseif kind == TreeSim.UnsampledUnary
                spec.allow_unsampled_unary ||
                    error("UnsampledUnary node at index $i is not valid in a reconstructed tree.")
                infectee = tree.host[tree.left[i]]
                _consume!(transmission_counts, (infectee, host, time)) ||
                    error("UnsampledUnary node at index $i does not match Transmission event.")
            end

            get(event_counts, (EpiSim.EK_Removal, host, time), 0) > 0 &&
                error("Removal event appears in tree (host=$host, time=$time).")

            _has_sampled_descendant(tree, i) ||
                error("Node at index $i has no sampled descendant.")

            key = (kind, host, time)
            key in seen && error("Duplicate node-event mapping detected for host=$host, time=$time.")
            push!(seen, key)
        end
    end

    return true
end

function _validate_bridge_log_contract(log::EventLog)
    sampled_signatures = Set{Tuple{EpiSim.EventKind,Int,Float64}}()

    for i in 1:length(log)
        kind = log.kind[i]
        _is_sampling(kind) || continue

        key = (kind, log.host[i], log.time[i])
        key in sampled_signatures &&
            error("Ambiguous duplicate retained sampling event signature: kind=$kind, host=$(log.host[i]), time=$(log.time[i]). EventLog-to-Tree conversion requires unique retained sampling signatures because Tree nodes do not store source event ids.")
        push!(sampled_signatures, key)
    end

    return nothing
end

function _validate_single_tree_contract(tree::TreeSim.Tree)
    isempty(tree) && return nothing
    _validate_component_contract(tree)
    return nothing
end

function _validate_component_contract(tree::TreeSim.Tree)
    root_ids = TreeSim.roots(tree)
    length(root_ids) == 1 ||
        error("Expected one rooted sampled-ancestry component per Tree; found $(length(root_ids)) roots. Use forest_from_eventlog for multi-component EventLogs.")

    for parent in eachindex(tree)
        for child in TreeSim.children(tree, parent)
            tree.time[parent] < tree.time[child] ||
                error("Equal-time or backward retained ancestry is not representable as a canonical TreeSim edge: edge $parent -> $child has parent time $(tree.time[parent]) and child time $(tree.time[child]). Retained ancestry times must be strictly increasing.")
        end
    end

    return nothing
end

function _increment!(counts::Dict{K,Int}, key::K) where {K}
    counts[key] = get(counts, key, 0) + 1
    return counts
end

function _consume!(counts::Dict{K,Int}, key::K) where {K}
    remaining = get(counts, key, 0)
    remaining == 0 && return false
    remaining == 1 ? delete!(counts, key) : (counts[key] = remaining - 1)
    return true
end

function _consume_any!(counts::Dict{K,Int}, keys) where {K}
    for key in keys
        _consume!(counts, key) && return true
    end
    return false
end

function _extract_sampled_forest(log::EventLog, spec::ExtractionSpec)
    tree = _extract_sampled_tree(log, spec)
    isempty(tree) && return TreeSim.Tree[]

    root_ids = sort(TreeSim.roots(tree); by=i -> (tree.time[i], tree.host[i], i))
    forest = [_subtree_as_tree(tree, root_id) for root_id in root_ids]

    for component in forest
        _validate_component_contract(component)
    end

    return forest
end

function _extract_sampled_tree(log::EventLog, spec::ExtractionSpec)
    if length(log) == 0
        return TreeSim.Tree()
    end

    max_host = max(maximum(log.host), maximum(log.infector))
    active = zeros(Int, max_host)
    tree = TreeSim.Tree()

    for i in length(log):-1:1
        kind = log.kind[i]
        time = log.time[i]
        host = log.host[i]

        if kind == EpiSim.EK_SerialSampling
            active[host] == 0 || error("SerialSampling encountered but host $host already has an active sampled lineage.")
            active[host] = _push_node!(tree, time, host, 0, 0, TreeSim.SampledLeaf)
        elseif kind == EpiSim.EK_FossilizedSampling
            child = active[host]
            if child == 0
                active[host] = _push_node!(tree, time, host, 0, 0, TreeSim.SampledLeaf)
            else
                node = _push_node!(tree, time, host, child, 0, TreeSim.SampledUnary)
                tree.parent[child] = node
                active[host] = node
            end
        elseif kind == EpiSim.EK_Transmission
            infectee = host
            infector = log.infector[i]
            child_node = active[infectee]
            child_node == 0 && continue

            parent_node = active[infector]
            if parent_node == 0
                node = _push_node!(tree, time, infector, child_node, 0, TreeSim.UnsampledUnary)
                tree.parent[child_node] = node
                active[infector] = node
            else
                node = _push_node!(tree, time, infector, parent_node, child_node, TreeSim.Binary)
                tree.parent[parent_node] = node
                tree.parent[child_node] = node
                active[infector] = node
            end
            active[infectee] = 0
        elseif kind == EpiSim.EK_Seeding
            child = active[host]
            child == 0 && continue
            node = _push_node!(tree, time, host, child, 0, TreeSim.Root)
            tree.parent[child] = node
            active[host] = 0
        end
    end

    any(!iszero, active) && error("Unresolved retained sampled lineages remain after processing the EventLog. The log is missing seeding ancestry for at least one sampled component.")
    isempty(tree) && return tree
    tree = _canonicalize(tree)
    spec.collapse_unsampled_unary && (tree = _collapse_unsampled_unary(tree))
    return tree
end

function _collapse_unsampled_unary(tree::TreeSim.Tree)
    drop = tree.kind .== TreeSim.UnsampledUnary
    any(drop) || return tree

    old_to_new = zeros(Int, length(tree))
    next = 0
    for i in eachindex(tree)
        drop[i] && continue
        next += 1
        old_to_new[i] = next
    end

    left = Int[]
    right = Int[]
    parent = Int[]

    for old in eachindex(tree)
        drop[old] && continue
        children = [_nearest_retained_descendant(tree, child, drop) for child in TreeSim.children(tree, old)]
        filter!(!iszero, children)
        push!(left, isempty(children) ? 0 : old_to_new[children[1]])
        push!(right, length(children) < 2 ? 0 : old_to_new[children[2]])
        retained_parent = _nearest_retained_parent(tree, old, drop)
        push!(parent, retained_parent == 0 ? 0 : old_to_new[retained_parent])
    end

    kept = findall(!, drop)
    return TreeSim.Tree(
        tree.time[kept],
        left,
        right,
        parent,
        tree.kind[kept],
        tree.host[kept],
        tree.label[kept],
    )
end

function _nearest_retained_descendant(tree::TreeSim.Tree, node::Int, drop::AbstractVector{Bool})
    node == 0 && return 0
    while drop[node]
        children = TreeSim.children(tree, node)
        isempty(children) && return 0
        node = only(children)
    end
    return node
end

function _nearest_retained_parent(tree::TreeSim.Tree, node::Int, drop::AbstractVector{Bool})
    parent = tree.parent[node]
    while parent != 0 && drop[parent]
        parent = tree.parent[parent]
    end
    return parent
end

function _subtree_as_tree(tree::TreeSim.Tree, root_id::Int)
    old_nodes = sort(collect(TreeSim.preorder(tree, root_id)))
    old_to_new = Dict{Int,Int}()
    for (new, old) in pairs(old_nodes)
        old_to_new[old] = new
    end

    return TreeSim.Tree(
        tree.time[old_nodes],
        [_remap_component(tree.left[old], old_to_new) for old in old_nodes],
        [_remap_component(tree.right[old], old_to_new) for old in old_nodes],
        [_remap_component(tree.parent[old], old_to_new) for old in old_nodes],
        tree.kind[old_nodes],
        tree.host[old_nodes],
        tree.label[old_nodes],
    )
end

function _push_node!(tree::TreeSim.Tree, time::Float64, host::Int, left::Int, right::Int, kind::TreeSim.NodeKind)
    push!(tree.time, time)
    push!(tree.left, left)
    push!(tree.right, right)
    push!(tree.parent, 0)
    push!(tree.kind, kind)
    push!(tree.host, host)
    push!(tree.label, kind == TreeSim.SampledLeaf ? host : 0)
    return length(tree)
end

function _canonicalize(tree::TreeSim.Tree)
    order = sortperm(eachindex(tree); by=i -> (tree.time[i], _kind_order(tree.kind[i]), i))
    old_to_new = zeros(Int, length(tree))
    for (new, old) in pairs(order)
        old_to_new[old] = new
    end

    return TreeSim.Tree(
        tree.time[order],
        [_remap(tree.left[old], old_to_new) for old in order],
        [_remap(tree.right[old], old_to_new) for old in order],
        [_remap(tree.parent[old], old_to_new) for old in order],
        tree.kind[order],
        tree.host[order],
        tree.label[order],
    )
end

_remap(i::Int, old_to_new::Vector{Int}) = i == 0 ? 0 : old_to_new[i]
_remap_component(i::Int, old_to_new::Dict{Int,Int}) = get(old_to_new, i, 0)
_is_sampling(kind::EpiSim.EventKind) = kind == EpiSim.EK_SerialSampling || kind == EpiSim.EK_FossilizedSampling
_kind_order(kind::TreeSim.NodeKind) = kind == TreeSim.Root ? 1 : kind == TreeSim.Binary || kind == TreeSim.UnsampledUnary || kind == TreeSim.SampledUnary ? 2 : 3

function _has_sampled_descendant(tree::TreeSim.Tree, i::Int)
    stack = [i]
    while !isempty(stack)
        j = pop!(stack)
        tree.kind[j] == TreeSim.SampledLeaf && return true
        tree.left[j] != 0 && push!(stack, tree.left[j])
        tree.right[j] != 0 && push!(stack, tree.right[j])
    end
    return false
end

end
