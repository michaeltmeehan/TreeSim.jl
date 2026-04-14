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

@testset "TreeSim canonical core" begin
    @testset "valid canonical binary tree" begin
        tree = canonical_binary_tree()

        @test validate_tree(tree)
        @test children(tree, 1) == [2, 3]
        @test children(tree, 2) == [4, 5]
        @test parent(tree, 4) == 2
        @test isbinary(tree, 1)
        @test isbinary(tree, 2)
        @test !isunary(tree, 2)
        @test isleaf(tree, 3)
        @test isinternal(tree, 1)
        @test leaves(tree) == [3, 4, 5]
        @test ancestors(tree, 5) == [2, 1]
        @test descendants(tree, 1) == [2, 4, 5, 3]
        @test descendants(tree, 2) == [4, 5]
    end

    @testset "valid canonical unary tree" begin
        tree = canonical_unary_tree()

        @test validate_tree(tree)
        @test isunary(tree, 1)
        @test isunary(tree, 2)
        @test isunary(tree, 3)
        @test leaves(tree) == [4]
        @test ancestors(tree, 4) == [3, 2, 1]
        @test descendants(tree, 1) == [2, 3, 4]
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
