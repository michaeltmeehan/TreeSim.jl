"""
    isleaf(node::DataFrameRow) -> Bool

Check if a node is a leaf node in a phylogenetic tree.

# Arguments
- `node::DataFrameRow`: A row from a DataFrame representing a node in the phylogenetic tree. The row should contain columns `left` and `right` indicating the left and right child nodes.

# Returns
- `Bool`: `true` if the node is a leaf node (i.e., both `left` and `right` are 0), `false` otherwise.
"""
function Phylo.isleaf(node::DataFrameRow)::Bool
    return node.left == 0 && node.right == 0
end


"""
    isbinary(node::DataFrameRow) -> Bool

Check if a node is a binary node in a phylogenetic tree.

# Arguments
- `node::DataFrameRow`: A row from a DataFrame representing a node in the phylogenetic tree. The row should contain a column `right` indicating the right child node.

# Returns
- `Bool`: `true` if the node is a binary node (i.e., `right` is not 0), `false` otherwise.
"""
function isbinary(node::DataFrameRow)::Bool
    return node.right != 0
end


"""
    isunary(node::DataFrameRow) -> Bool

Check if a node is a unary node in a phylogenetic tree.

# Arguments
- `node::DataFrameRow`: A row from a DataFrame representing a node in the phylogenetic tree. The row should contain columns `left` and `right` indicating the left and right child nodes.

# Returns
- `Bool`: `true` if the node is a unary node (i.e., `left` is not 0 and `right` is 0), `false` otherwise.
"""
function isunary(node::DataFrameRow)::Bool
    return node.left != 0 && node.right == 0
end


"""
    isroot(node::DataFrameRow) -> Bool

Check if a node is the root node in a phylogenetic tree.

# Arguments
- `node::DataFrameRow`: A row from a DataFrame representing a node in the phylogenetic tree. The row should contain a column `host` indicating the ID of the host node.

# Returns
- `Bool`: `true` if the node is the root node (i.e., `host` is 0), `false` otherwise.
"""
function Phylo.isroot(node::DataFrameRow)::Bool
    return node.host == 0
end


"""
    binarize_tree(tree::DataFrame) -> DataFrame

Transform a phylogenetic tree into its binary form by collapsing unary nodes and retaining only leaf, binary, and root nodes.

# Arguments
- `tree::DataFrame`: A DataFrame representing the phylogenetic tree. The DataFrame should contain columns `id`, `left`, `right`, and `host`.

# Returns
- `DataFrame`: A new DataFrame representing the binary form of the input phylogenetic tree. The new DataFrame will have updated `id`, `left`, and `right` columns to reflect the collapsed unary nodes.
"""
function binarize_tree(tree::DataFrame)
    ctree = copy(tree)
    binary_tree = DataFrame()
    for row in eachrow(ctree)
        if isleaf(row)
            push!(binary_tree, row)
        elseif isbinary(row)
            left = row.left
            right = row.right
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            while isunary(tree[right, :])
                right = tree[right, :left]
            end
            row.left = left
            row.right = right
            push!(binary_tree, row)
        elseif isroot(row)
            left = row.left
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            row.left = left
            push!(binary_tree, row)
        end
    end
    new_ids = Dict(zip(binary_tree.id, 1:nrow(binary_tree)))
    new_ids[0] = 0
    binary_tree.id = [new_ids[id] for id in binary_tree.id]
    binary_tree.left = [new_ids[id] for id in binary_tree.left]
    binary_tree.right = [new_ids[id] for id in binary_tree.right]
    return binary_tree
end


"""
    relabel!(wtree::DataFrame) -> DataFrame

Assign unique labels to each node in a phylogenetic tree using Cantor pairing.

# Arguments
- `wtree::DataFrame`: A DataFrame representing the phylogenetic tree. The DataFrame should contain columns `id`, `left`, `right`, `leaf_id`, and `host`.

# Returns
- `DataFrame`: The modified DataFrame with updated node labels.
"""
function relabel!(wtree::DataFrame)::DataFrame
    new_ids = Dict{Int64, Int64}()
    for node in eachrow(wtree)
        if node.left == 0 && node.right == 0 # Leaf node
            if node.leaf_id == node.host # Sampled leaf
                new_ids[node.id] = CantorPair(node.id, node.host)
                node.id = new_ids[node.id]
            else # Transmission leaf
                new_ids[node.id] = CantorPair(0, node.leaf_id)
                node.id = new_ids[node.id]
            end
        else    # Internal node
            new_ids[node.id] = CantorPair(node.id, node.host)
            node.id = new_ids[node.id]
            node.left = new_ids[node.left]
            node.right = node.right == 0 ? 0 : new_ids[node.right]
        end
    end
    return wtree
end


"""
    normalize!(tree::DataFrame) -> DataFrame

Sequentially label nodes in a phylogenetic tree.

# Arguments
- `tree::DataFrame`: A DataFrame representing the phylogenetic tree. The DataFrame should contain columns `id`, `left`, and `right`.

# Returns
- `DataFrame`: The modified DataFrame with nodes labeled sequentially.
"""
function normalize!(tree::DataFrame)::DataFrame
    new_labels = Dict(zip(tree.id, 1:nrow(tree)))
    new_labels[0] = 0
    tree.id = 1:nrow(tree)
    tree.left = [new_labels[node] for node in tree.left]
    tree.right = [new_labels[node] for node in tree.right]
    return tree
end


"""
    join_wtrees(trees::Vector{DataFrame}) -> DataFrame

Combine multiple within-host phylogenetic trees into a single tree, remove duplicate nodes, and normalize node labels.

# Arguments
- `trees::Vector{DataFrame}`: A vector of DataFrames, each representing a within-host phylogenetic tree.

# Returns
- `DataFrame`: A single DataFrame representing the combined phylogenetic tree with duplicate nodes removed and node labels normalized.
"""
function join_wtrees(trees::Vector{DataFrame})::DataFrame
    tree = vcat(trees...)
    tree = tree[(tree.leaf_id .== 0 .|| tree.leaf_id .== tree.host), :] # Remove duplicate nodes
    sort!(tree, :t, rev=true)
    normalize!(tree)
    return tree
end