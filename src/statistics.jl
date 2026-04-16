"""
    TreeStats

Summary statistics for a sampled-ancestry `Tree`.
"""
struct TreeStats
    ntips::Int
    height::Float64
    total_branch_length::Float64
    mean_internal_branch::Float64
    mean_terminal_branch::Float64
    internal_terminal_ratio::Float64
    mean_node_depth::Float64
    var_node_depth::Float64
    sackin::Int
    colless::Int
    cherries::Int
    branch_length_cv::Float64
    subtree_size_cv::Float64
    mean_child_imbalance::Float64
    ladder_fraction::Float64
    root_to_tip_cv::Float64
end

"""
    tree_statistics(tree::Tree)

Compute legacy outbreak tree summary statistics on a canonical TreeSim tree.

Unary chains are collapsed for branch-length summaries, matching the supplied
legacy workflow. Topological imbalance summaries are computed over binary
nodes. Empty trees and trees without sampled tips are rejected.
"""
function tree_statistics(tree::Tree)
    isempty(tree) && throw(ArgumentError("tree must be non-empty."))
    validate_tree(tree; require_single_root=true, require_reachable=true)

    n = length(tree)
    subtree_size = zeros(Int, n)
    edge_depth = zeros(Int, n)

    total_branch_length = 0.0
    internal_sum = 0.0
    terminal_sum = 0.0
    internal_count = 0
    terminal_count = 0
    branch_sum = 0.0
    branch_sqsum = 0.0
    branch_count = 0

    ntips = 0
    root_time = Inf
    max_time = -Inf
    colless = 0
    cherries = 0
    subtree_sum = 0.0
    subtree_sqsum = 0.0
    subtree_count = 0
    imbalance_sum = 0.0
    imbalance_count = 0
    ladder_nodes = 0
    binary_nodes = 0

    for i in postorder(tree)
        t = tree.time[i]
        max_time = max(max_time, t)
        tree.parent[i] == 0 && (root_time = min(root_time, t))

        if tree.kind[i] == SampledLeaf
            subtree_size[i] = 1
            ntips += 1
        end

        l = tree.left[i]
        r = tree.right[i]
        l != 0 && (subtree_size[i] += subtree_size[l])
        r != 0 && (subtree_size[i] += subtree_size[r])

        if tree.kind[i] == Binary
            sl = subtree_size[l]
            sr = subtree_size[r]
            colless += abs(sl - sr)
            (sl == 1 && sr == 1) && (cherries += 1)

            s = subtree_size[i]
            subtree_sum += s
            subtree_sqsum += s * s
            subtree_count += 1

            imbalance_sum += abs(sl - sr) / (sl + sr)
            imbalance_count += 1
            min(sl, sr) == 1 && (ladder_nodes += 1)
            binary_nodes += 1
        end

        if tree.parent[i] != 0 && tree.kind[i] != UnsampledUnary && tree.kind[i] != SampledUnary
            p = tree.parent[i]
            bl = tree.time[i] - tree.time[p]

            while p != 0 && (tree.kind[p] == UnsampledUnary || tree.kind[p] == SampledUnary)
                gp = tree.parent[p]
                gp == 0 && break
                bl += tree.time[p] - tree.time[gp]
                p = gp
            end

            total_branch_length += bl
            branch_sum += bl
            branch_sqsum += bl * bl
            branch_count += 1

            if tree.kind[i] == SampledLeaf
                terminal_sum += bl
                terminal_count += 1
            elseif tree.kind[i] == Binary
                internal_sum += bl
                internal_count += 1
            end
        end
    end

    ntips > 0 || throw(ArgumentError("tree must contain at least one sampled tip."))

    sackin = 0
    for i in preorder(tree)
        for child in children(tree, i)
            edge_depth[child] = edge_depth[i] + 1
        end
        tree.kind[i] == SampledLeaf && (sackin += edge_depth[i])
    end

    depth_sum = 0.0
    depth_sqsum = 0.0
    depth_count = 0
    for i in eachindex(tree)
        if tree.kind[i] == Binary
            d = max_time - tree.time[i]
            depth_sum += d
            depth_sqsum += d * d
            depth_count += 1
        end
    end

    mean_node_depth = _safe_mean(depth_sum, depth_count)
    var_node_depth = depth_count == 0 ? NaN : max(0.0, depth_sqsum / depth_count - mean_node_depth^2)
    mean_internal = _safe_mean(internal_sum, internal_count)
    mean_terminal = _safe_mean(terminal_sum, terminal_count)
    branch_mean = _safe_mean(branch_sum, branch_count)
    branch_var = branch_count == 0 ? NaN : max(0.0, branch_sqsum / branch_count - branch_mean^2)
    mean_subtree = _safe_mean(subtree_sum, subtree_count)
    var_subtree = subtree_count == 0 ? NaN : max(0.0, subtree_sqsum / subtree_count - mean_subtree^2)

    rtt_sum = 0.0
    rtt_sqsum = 0.0
    for i in eachindex(tree)
        if tree.kind[i] == SampledLeaf
            d = tree.time[i] - root_time
            rtt_sum += d
            rtt_sqsum += d * d
        end
    end
    mean_rtt = rtt_sum / ntips
    var_rtt = max(0.0, rtt_sqsum / ntips - mean_rtt^2)

    return TreeStats(
        ntips,
        max_time - root_time,
        total_branch_length,
        mean_internal,
        mean_terminal,
        _safe_ratio(mean_internal, mean_terminal),
        mean_node_depth,
        var_node_depth,
        sackin,
        colless,
        cherries,
        _safe_ratio(sqrt(branch_var), branch_mean),
        _safe_ratio(sqrt(var_subtree), mean_subtree),
        _safe_mean(imbalance_sum, imbalance_count),
        _safe_ratio(ladder_nodes, binary_nodes),
        _safe_ratio(sqrt(var_rtt), mean_rtt),
    )
end

_safe_mean(x, n::Integer) = n == 0 ? NaN : x / n
_safe_ratio(x, y) = iszero(y) || isnan(y) ? NaN : x / y
