using Test
using TreeSim

function canonical_binary_tree()
    return Tree(
        [0.0, 0.4, 0.8, 1.1, 1.3],
        [2, 4, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [0, 1, 1, 2, 2],
        [Root, Binary, SampledLeaf, SampledLeaf, SampledLeaf],
        [0, 0, 0, 0, 0],
        [0, 0, 101, 102, 103],
    )
end

function canonical_unary_tree()
    return Tree(
        [0.0, 0.4, 0.9, 1.5],
        [2, 3, 4, 0],
        [0, 0, 0, 0],
        [0, 1, 2, 3],
        [Root, UnsampledUnary, SampledUnary, SampledLeaf],
        [1, 2, 2, 2],
        [0, 0, 301, 302],
    )
end

function two_cherry_tree()
    return Tree(
        [0.0, 0.3, 0.4, 0.8, 0.9, 1.0, 1.1],
        [2, 4, 6, 0, 0, 0, 0],
        [3, 5, 7, 0, 0, 0, 0],
        [0, 1, 1, 2, 2, 3, 3],
        [Root, Binary, Binary, SampledLeaf, SampledLeaf, SampledLeaf, SampledLeaf],
        [0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 401, 402, 403, 404],
    )
end

@testset "TreeSim canonical core" begin
    @testset "valid canonical binary tree" begin
        tree = canonical_binary_tree()

        @test validate_tree(tree)
        @test roots(tree) == [1]
        @test root(tree) == 1
        @test children(tree, 1) == [2, 3]
        @test children(tree, 2) == [4, 5]
        @test parent(tree, 4) == 2
        @test isbinary(tree, 1)
        @test isbinary(tree, 2)
        @test !isunary(tree, 2)
        @test isleaf(tree, 3)
        @test isleaf(tree, 4)
        @test !isleaf(tree, 1)
        @test !isleaf(tree, 2)
        @test isinternal(tree, 1)
        @test leaves(tree) == [3, 4, 5]
        @test tips(tree) == [4, 5, 3]
        @test length(tips(tree)) == nleaves(tree)
        @test internal_nodes(tree) == [1, 2]
        @test nnodes(tree) == 5
        @test nleaves(tree) == 3
        @test ninternal(tree) == 2
        @test ancestors(tree, 5) == [2, 1]
        @test descendants(tree, 1) == [2, 4, 5, 3]
        @test descendants(tree, 2) == [4, 5]
        @test subtree_nodes(tree, 3) == [3]
        @test subtree_nodes(tree, 2) == [2, 4, 5]
        @test subtree_nodes(tree, 1) == [1, 2, 4, 5, 3]
        @test node_depths(tree) ≈ [0.0, 0.4, 0.8, 1.1, 1.3]
        @test root_to_tip_distances(tree) ≈ [1.1, 1.3, 0.8]
        @test tree_height(tree) ≈ 1.3
        @test mean_root_to_tip_distance(tree) ≈ sum(root_to_tip_distances(tree)) / nleaves(tree)
        @test ncherries(tree) == 1
    end

    @testset "valid canonical unary tree" begin
        tree = canonical_unary_tree()

        @test validate_tree(tree)
        @test root(tree) == 1
        @test isunary(tree, 1)
        @test isunary(tree, 2)
        @test isunary(tree, 3)
        @test !isleaf(tree, 1)
        @test isleaf(tree, 4)
        @test leaves(tree) == [4]
        @test tips(tree) == [4]
        @test length(tips(tree)) == nleaves(tree)
        @test internal_nodes(tree) == [1, 2, 3]
        @test nnodes(tree) == 4
        @test nleaves(tree) == 1
        @test ninternal(tree) == 3
        @test ancestors(tree, 4) == [3, 2, 1]
        @test descendants(tree, 1) == [2, 3, 4]
        @test subtree_nodes(tree, 4) == [4]
        @test subtree_nodes(tree, 2) == [2, 3, 4]
        @test subtree_nodes(tree, 1) == [1, 2, 3, 4]
        @test node_depths(tree) ≈ [0.0, 0.4, 0.9, 1.5]
        @test root_to_tip_distances(tree) ≈ [1.5]
        @test tree_height(tree) ≈ 1.5
        @test mean_root_to_tip_distance(tree) ≈ 1.5
        @test ncherries(tree) == 0
    end

    @testset "depth-oriented summaries" begin
        tree = two_cherry_tree()

        @test validate_tree(tree)
        @test tips(tree) == [4, 5, 6, 7]
        @test node_depths(tree) ≈ [0.0, 0.3, 0.4, 0.8, 0.9, 1.0, 1.1]
        @test length(node_depths(tree)) == nnodes(tree)
        @test node_depths(tree)[tips(tree)] ≈ root_to_tip_distances(tree)
        @test root_to_tip_distances(tree) ≈ [0.8, 0.9, 1.0, 1.1]
        @test tree_height(tree) ≈ maximum(root_to_tip_distances(tree))
        @test mean_root_to_tip_distance(tree) ≈ 0.95
        @test ncherries(tree) == 2
    end

    @testset "traversal iterators" begin
        tree = canonical_binary_tree()

        @test collect(preorder(tree)) == [1, 2, 4, 5, 3]
        @test collect(preorder(tree, 2)) == [2, 4, 5]
        @test collect(postorder(tree)) == [4, 5, 2, 3, 1]
        @test collect(postorder(tree, 2)) == [4, 5, 2]
        @test collect(breadthfirst(tree)) == [1, 2, 3, 4, 5]
        @test collect(breadthfirst(tree, 2)) == [2, 4, 5]

        unary = canonical_unary_tree()
        @test collect(preorder(unary)) == [1, 2, 3, 4]
        @test collect(postorder(unary)) == [4, 3, 2, 1]
        @test collect(breadthfirst(unary)) == [1, 2, 3, 4]
    end

    @testset "root discovery handles provisional root counts" begin
        empty = Tree()
        @test roots(empty) == Int[]
        @test_throws ErrorException root(empty)
        @test collect(preorder(empty)) == Int[]
        @test collect(postorder(empty)) == Int[]
        @test collect(breadthfirst(empty)) == Int[]
        @test tips(empty) == Int[]
        @test nnodes(empty) == 0
        @test nleaves(empty) == 0
        @test ninternal(empty) == 0
        @test node_depths(empty) == Float64[]
        @test root_to_tip_distances(empty) == Float64[]
        @test tree_height(empty) == 0.0
        @test mean_root_to_tip_distance(empty) == 0.0
        @test ncherries(empty) == 0

        multi_root_tree = Tree(
            [0.0, 0.1, 0.5, 0.8],
            [3, 4, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 1, 2],
            [Root, Root, SampledLeaf, SampledLeaf],
            [0, 0, 0, 0],
            [0, 0, 201, 202],
        )

        @test roots(multi_root_tree) == [1, 2]
        @test_throws ErrorException root(multi_root_tree)
        @test collect(preorder(multi_root_tree)) == [1, 3, 2, 4]
        @test collect(postorder(multi_root_tree)) == [3, 1, 4, 2]
        @test collect(breadthfirst(multi_root_tree)) == [1, 2, 3, 4]
        @test tips(multi_root_tree) == [3, 4]
        @test length(tips(multi_root_tree)) == nleaves(multi_root_tree)
        @test collect(preorder(multi_root_tree, 2)) == [2, 4]
        @test collect(postorder(multi_root_tree, 2)) == [4, 2]
        @test collect(breadthfirst(multi_root_tree, 2)) == [2, 4]
        @test subtree_nodes(multi_root_tree, 3) == [3]
        @test subtree_nodes(multi_root_tree, 1) == [1, 3]
        @test subtree_nodes(multi_root_tree, 2) == [2, 4]
        @test node_depths(multi_root_tree) ≈ [0.0, 0.0, 0.5, 0.7]
        @test node_depths(multi_root_tree)[tips(multi_root_tree)] ≈ root_to_tip_distances(multi_root_tree)
        @test root_to_tip_distances(multi_root_tree) ≈ [0.5, 0.7]
        @test tree_height(multi_root_tree) ≈ 0.7
        @test mean_root_to_tip_distance(multi_root_tree) ≈ 0.6
        @test ncherries(multi_root_tree) == 0
    end

    @testset "node row view and branch length" begin
        tree = canonical_binary_tree()
        node = tree[4]

        @test node isa Node
        @test node.id == 4
        @test node.time == 1.1
        @test node.parent == 2
        @test node.kind == SampledLeaf
        @test branch_length(tree, 2) == 0.4
        @test branch_length(tree, 4) ≈ 0.7
        @test_throws ErrorException branch_length(tree, 1)
    end

    @testset "canonical validation requires one reachable root by default" begin
        multi_root_tree = Tree(
            [0.0, 0.1, 0.5, 0.8],
            [3, 4, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 1, 2],
            [Root, Root, SampledLeaf, SampledLeaf],
            [0, 0, 0, 0],
            [0, 0, 201, 202],
        )

        @test_throws ErrorException validate_tree(multi_root_tree)
        @test validate_tree(multi_root_tree; require_single_root=false)
    end

    @testset "root identity and node-kind degrees are explicit" begin
        nonroot_without_parent = Tree(
            [0.0],
            [0],
            [0],
            [0],
            [SampledLeaf],
            [0],
            [0],
        )

        root_with_parent = Tree(
            [0.0, 1.0],
            [2, 0],
            [0, 0],
            [0, 1],
            [Root, Root],
            [0, 0],
            [0, 0],
        )

        binary_with_one_child = Tree(
            [0.0, 1.0],
            [2, 0],
            [0, 0],
            [0, 1],
            [Root, Binary],
            [0, 0],
            [0, 0],
        )

        leaf_with_child = Tree(
            [0.0, 0.5, 1.0],
            [2, 3, 0],
            [0, 0, 0],
            [0, 1, 2],
            [Root, SampledLeaf, SampledLeaf],
            [0, 0, 0],
            [0, 0, 0],
        )

        @test_throws ErrorException validate_tree(nonroot_without_parent)
        @test_throws ErrorException validate_tree(root_with_parent)
        @test_throws ErrorException validate_tree(binary_with_one_child)
        @test_throws ErrorException validate_tree(leaf_with_child)
    end

    @testset "parent-child references are reciprocal and unique" begin
        child_missing_parent = Tree(
            [0.0, 1.0],
            [2, 0],
            [0, 0],
            [0, 0],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        parent_missing_child = Tree(
            [0.0, 1.0],
            [0, 0],
            [0, 0],
            [0, 1],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        same_child_twice = Tree(
            [0.0, 1.0],
            [2, 0],
            [2, 0],
            [0, 1],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        shared_child = Tree(
            [0.0, 0.4, 0.8, 1.2],
            [2, 4, 0, 0],
            [3, 4, 0, 0],
            [0, 1, 1, 2],
            [Root, Binary, SampledLeaf, SampledLeaf],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        )

        @test_throws ErrorException validate_tree(child_missing_parent)
        @test_throws ErrorException validate_tree(parent_missing_child)
        @test_throws ErrorException validate_tree(same_child_twice)
        @test_throws ErrorException validate_tree(shared_child)
    end

    @testset "invalid indices and cycles are rejected" begin
        invalid_child_index = Tree(
            [0.0, 1.0],
            [3, 0],
            [0, 0],
            [0, 1],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        invalid_parent_index = Tree(
            [0.0, 1.0],
            [2, 0],
            [0, 0],
            [0, 3],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        self_child = Tree(
            [0.0],
            [1],
            [0],
            [0],
            [Root],
            [0],
            [0],
        )

        parent_cycle = Tree(
            [0.0, 1.0, 2.0],
            [2, 3, 0],
            [0, 0, 0],
            [0, 3, 2],
            [Root, UnsampledUnary, UnsampledUnary],
            [0, 0, 0],
            [0, 0, 0],
        )

        @test_throws ErrorException validate_tree(invalid_child_index)
        @test_throws ErrorException validate_tree(invalid_parent_index)
        @test_throws ErrorException validate_tree(self_child)
        @test_throws ErrorException validate_tree(parent_cycle; require_reachable=false)
    end

    @testset "time and storage ordering are canonical invariants" begin
        equal_edge_time = Tree(
            [0.0, 0.0],
            [2, 0],
            [0, 0],
            [0, 1],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )

        decreasing_storage_time = Tree(
            [0.0, 1.0, 0.5],
            [3, 0, 0],
            [2, 0, 0],
            [0, 1, 1],
            [Root, SampledLeaf, SampledLeaf],
            [0, 0, 0],
            [0, 0, 0],
        )

        backward_child_index = Tree(
            [0.0, 0.5, 1.0],
            [0, 1, 0],
            [3, 0, 0],
            [2, 0, 2],
            [SampledLeaf, Root, SampledLeaf],
            [0, 0, 0],
            [0, 0, 0],
        )

        @test_throws ErrorException validate_tree(equal_edge_time)
        @test_throws ErrorException validate_tree(decreasing_storage_time)
        @test_throws ErrorException validate_tree(backward_child_index)
    end
end
