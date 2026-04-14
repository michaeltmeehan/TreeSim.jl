"""
    sample_bounded_times(times::Vector{Float64},
                         leaves::Vector{Int64},
                         Nₑ::Float64,
                         bound::Float64,
                         n_sam::Int64) -> Tuple{Matrix{Float64}, Vector{Float64}}

Sample coalescence times using the bounded coalescent model and calculate the likelihood of each sample.

# Arguments
- `times::Vector{Float64}`: A vector of time points.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.
- `n_sam::Int64`: The number of samples to generate.

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}}`: A tuple containing:
  - `Matrix{Float64}`: A matrix of sampled coalescence times.
  - `Vector{Float64}`: A vector of likelihoods for each sample.
"""
function sample_bounded_times(times::Vector{Float64},
                              leaves::Vector{Int64},
                              Nₑ::Float64,
                              bound::Float64,
                              n_sam::Int64)::Tuple{Matrix{Float64}, Vector{Float64}}

    # Total the leaves
    total_leaves = sum(leaves)

    # Storage for samples
    sampled_times = Matrix{Float64}(undef, n_sam, (total_leaves - 1))

    # Storage for likelihood
    likelihood = ones(n_sam)

    forward_probs = forward_algorithm(times, leaves, Nₑ, bound)

    # Random numbers for sampling
    u = rand(n_sam, (total_leaves - 1))

    lineages = fill(1,length(times) + 1)
    const_lower = Vector{Float64}(undef, total_leaves - 1)
    const_upper = Vector{Float64}(undef, total_leaves - 1)
    const_lineages = Vector{Int64}(undef, total_leaves - 1)
    const_events = Vector{Int64}(undef, total_leaves - 1)

    for n in 1:n_sam
        likelihood[n] *= backward_sampler!(lineages, forward_probs, times, leaves, Nₑ, bound)

        likelihood[n] *= constrain_coalescences!(const_lower, const_upper, const_lineages, const_events,
                                                    lineages, times, leaves, Nₑ, bound)

        for c in 1:(total_leaves - 1)
            dl = convert(Float64, const_lineages[c])
            z = (Nₑ / (dl - 1.0)) * (1.0 - exp(((dl - 1.0) / Nₑ) *
                    (const_lower[c] - const_upper[c])))

            sampled_times[n, c] = const_upper[c] + (Nₑ / (dl - 1.0)) *
                                    log(1.0 - ((dl - 1.0) / Nₑ) * z * u[n, c])

            likelihood[n] *= (1.0 / z) * exp(((dl - 1.0) / Nₑ) *
                                (sampled_times[n, c] - const_upper[c]))
        end
    end
  return sampled_times, likelihood
end


"""
    sample_bounded_times(times::Vector{Float64},
                         leaves::Vector{Int64},
                         Nₑ::Float64,
                         bound::Float64) -> Tuple{Vector{Float64}, Float64}

Sample coalescence times using the bounded coalescent model and calculate the likelihood for a single sample.

# Arguments
- `times::Vector{Float64}`: A vector of time points.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.

# Returns
- `Tuple{Vector{Float64}, Float64}`: A tuple containing:
  - `Vector{Float64}`: A vector of sampled coalescence times.
  - `Float64`: The likelihood of the sampled coalescence times.
"""
function sample_bounded_times(times::Vector{Float64},
                              leaves::Vector{Int64},
                              Nₑ::Float64,
                              bound::Float64)::Tuple{Vector{Float64}, Float64}
    sampled_times, likelihood = sample_bounded_times(times, leaves, Nₑ, bound, 1)
    return vec(sampled_times), likelihood[1]
end


"""
    sample_topology(leaf_times::Vector{Float64},
                    leaves::Vector{Int64},
                    coalescence_times::Vector{Float64}) -> Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}

Sample the topology of the phylogenetic tree based on coalescence times.

# Arguments
- `leaf_times::Vector{Float64}`: A vector of leaf (sampling) times.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `coalescence_times::Vector{Float64}`: A vector of coalescence times.

# Returns
- `Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}`: A tuple containing:
  - `Matrix{Int64}`: A matrix representing the edges of the phylogenetic tree.
  - `Vector{Float64}`: A vector of edge lengths.
  - `Float64`: The likelihood of the sampled tree topology.
  - `DataFrame`: A dataframe containing node information, including times, IDs, and parent-child relationships.
"""
function sample_topology(leaf_times::Vector{Float64},
                         leaves::Vector{Int64},
                         coalescence_times::Vector{Float64})::Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}

    # Total the leaves
    total_leaves = sum(leaves)

    # Internal node ancestors
    edge = Matrix{Int64}(undef, 2*(total_leaves - 1), 2)
    edge_length = Vector{Float64}(undef, 2*(total_leaves - 1))
    # Node ancestors and times
    node_times = fill(0., 2 * total_leaves - 1)

    # Counters
    l_index = total_leaves
    n_index = 2 * total_leaves - 1

    # Current information
    active_nodes = fill(0, 2*total_leaves - 1)
    total_active_nodes = 0

    # Node dataframe
    node_df = DataFrame(t = reduce(vcat, fill.(leaf_times, leaves)),
                        id = 1:total_leaves,
                        left = fill(0, total_leaves),
                        right = fill(0, total_leaves))

    # Random numbers for sampling
    u = rand(2*(total_leaves - 1))

    # Likelihood
    likelihood = 1.0

    k = length(leaf_times)

    for i in 1:leaves[k]
        active_nodes[l_index] = 1
        total_active_nodes += 1
        node_times[l_index] = leaf_times[k]
        l_index -= 1
    end

    k -= 1

    c = length(coalescence_times)

    while (c >= 1)

        if (k < 1 || leaf_times[k] < coalescence_times[c])
            anc_1 = 1
            sum_prob = convert(Float64, active_nodes[anc_1]) / convert(Float64, total_active_nodes)
            while (sum_prob < u[2 * c])
                anc_1 += 1
                sum_prob += convert(Float64, active_nodes[anc_1]) / convert(Float64, total_active_nodes)
            end

            edge[2 * c, 1] = n_index
            edge[2 * c, 2] = anc_1

            edge_length[2 * c] = node_times[anc_1] - coalescence_times[c]

            likelihood *= 2.0 / convert(Float64, total_active_nodes)

            active_nodes[anc_1] = 0
            total_active_nodes -= 1

            anc_2 = 1
            sum_prob = convert(Float64, active_nodes[anc_2]) / convert(Float64, total_active_nodes)
            while (sum_prob < u[2 * c - 1 ])
                anc_2 += 1
                sum_prob += convert(Float64, active_nodes[anc_2]) / convert(Float64, total_active_nodes)
            end
            edge[2 * c - 1, 1] = n_index
            edge[2 * c - 1, 2] = anc_2

            edge_length[2 * c - 1] = node_times[anc_2] - coalescence_times[c]

            likelihood *= 1.0 / convert(Float64, total_active_nodes)

            active_nodes[anc_2] = 0

            active_nodes[n_index] = 1
            node_times[n_index] = coalescence_times[c]

            push!(node_df, [node_times[n_index] n_index anc_1 anc_2])

            c -= 1
            n_index -= 1
        else
            for i in 1:leaves[k]
                active_nodes[l_index] = 1
                total_active_nodes += 1

                node_times[l_index] = leaf_times[k]
                l_index -= 1
            end
            k -= 1
        end
    end
    return edge, edge_length, likelihood, node_df
end


"""
    sample_wtree(host::Host, Nₑ::Float64; relabel=true) -> DataFrame

Generate a within-host phylogenetic tree for a given host based on their leaf (sampling) times and other characteristics.

# Arguments
- `host::Host`: A `Host` object containing information about the host, including leaf times, root time, and IDs.
- `Nₑ::Float64`: The effective population size.
- `relabel::Bool`: A boolean flag indicating whether to relabel the nodes in the resulting tree (default is true).

# Returns
- `DataFrame`: A dataframe representing the within-host phylogenetic tree, containing columns for times (`t`), node IDs (`id`), left and right child IDs (`left`, `right`), leaf IDs (`leaf_id`), host IDs (`host`), and node types (`type`).
"""
function sample_wtree(host::Host, Nₑ::Float64; relabel=true)
    isnothing(host.infector) && return DataFrame(t = [], id = [], left = [], right = [], leaf_id = [], host = [], type = [])
    leaves = fill(1, length(host.leaf_times))
    tot_leaves = length(leaves)
    if tot_leaves == 1 
        wtree = DataFrame(t = [host.leaf_times[1], host.root], id = [1, 0], left = [0, 1], right = [0, 0])
    else
        _, _, _, wtree = sample_wtree(host.leaf_times, leaves, Nₑ, host.root)
        push!(wtree, (host.root, 0, tot_leaves + 1, 0))
    end
    wtree.leaf_id .= vcat(host.leaf_ids, fill(0, tot_leaves))
    wtree.host .= host.id
    wtree.type .= host.type
    relabel && relabel!(wtree)
    wtree[end, :host] = host.infector
    wtree[end, :type] = host.infector_type
    return wtree
end


"""
    sample_wtree(leaf_times::Vector{Float64},
                 leaves::Vector{Int64},
                 Nₑ::Float64,
                 bound::Float64) -> Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}

Generate a phylogenetic tree based on leaf (sampling) times, number of leaves, effective population size, and a time bound.

# Arguments
- `leaf_times::Vector{Float64}`: A vector of leaf (sampling) times.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.

# Returns
- `Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}`: A tuple containing:
  - `Matrix{Int64}`: A matrix representing the edges of the phylogenetic tree.
  - `Vector{Float64}`: A vector of edge lengths.
  - `Float64`: The likelihood of the sampled tree.
  - `DataFrame`: A dataframe containing node information, including times, IDs, and parent-child relationships.
"""
function sample_wtree(leaf_times::Vector{Float64},
                      leaves::Vector{Int64},
                      Nₑ::Float64,
                      bound::Float64)::Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}
    coalescence_times, time_likelihood = sample_bounded_times(leaf_times, leaves, Nₑ, bound)
    edge, edge_length, top_likelihood, node_df = sample_topology(leaf_times, leaves, coalescence_times)
    return edge, edge_length, time_likelihood * top_likelihood, node_df
end