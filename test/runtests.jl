using Test
using TreeSim

function binary_tree()
    return Tree(
        [0.0, 1.0, 1.2],
        [2, 0, 0],
        [3, 0, 0],
        [0, 1, 1],
        [Root, SampledLeaf, SampledLeaf],
        [0, 0, 0],
        [0, 101, 102],
    )
end

function unary_tree()
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

function multi_root_tree()
    return Tree(
        [0.0, 0.5, 0.1, 0.8],
        [2, 0, 4, 0],
        [0, 0, 0, 0],
        [0, 1, 0, 3],
        [Root, SampledLeaf, Root, SampledLeaf],
        [0, 0, 0, 0],
        [0, 201, 0, 202],
    )
end

@testset "TreeSim core" begin
    @testset "valid binary tree" begin
        tree = binary_tree()
        @test validate_tree(tree)
        @test validate_tree(tree; require_single_root=true)
        @test validate_tree(tree; require_reachable=true)
        @test children(tree, 1) == [2, 3]
        @test parent(tree, 2) == 1
        @test isbinary(tree, 1)
        @test !isunary(tree, 1)
        @test isleaf(tree, 2)
        @test isinternal(tree, 1)
        @test leaves(tree) == [2, 3]
        @test ancestors(tree, 2) == [1]
        @test descendants(tree, 1) == [2, 3]
    end

    @testset "valid unary tree" begin
        tree = unary_tree()
        @test validate_tree(tree)
        @test validate_tree(tree; require_single_root=true, require_reachable=true)
        @test isunary(tree, 1)
        @test isunary(tree, 2)
        @test isunary(tree, 3)
        @test leaves(tree) == [4]
        @test ancestors(tree, 4) == [3, 2, 1]
        @test descendants(tree, 1) == [2, 3, 4]
    end

    @testset "multi-root validation policy" begin
        tree = multi_root_tree()
        @test validate_tree(tree)
        @test validate_tree(tree; require_reachable=true)
        @test_throws ErrorException validate_tree(tree; require_single_root=true)
    end

    @testset "strict reachability catches orphaned nodes" begin
        tree = Tree(
            [0.0, 1.0, 1.1],
            [2, 0, 0],
            [0, 0, 0],
            [0, 1, 1],
            [Root, SampledLeaf, SampledLeaf],
            [0, 0, 0],
            [0, 101, 102],
        )
        @test validate_tree(tree)
        @test_throws ErrorException validate_tree(tree; require_reachable=true)
    end

    @testset "root identity is explicit" begin
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
        @test_throws ErrorException validate_tree(nonroot_without_parent)
        @test_throws ErrorException validate_tree(root_with_parent)
    end

    @testset "invalid parent-child consistency" begin
        tree = Tree(
            [0.0, 1.0],
            [2, 0],
            [0, 0],
            [0, 0],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )
        @test_throws ErrorException validate_tree(tree)
    end

    @testset "invalid time ordering" begin
        tree = Tree(
            [1.0, 0.5],
            [2, 0],
            [0, 0],
            [0, 1],
            [Root, SampledLeaf],
            [0, 0],
            [0, 0],
        )
        @test_throws ErrorException validate_tree(tree)
    end

    @testset "node row view and branch length" begin
        tree = binary_tree()
        node = tree[2]
        @test node isa Node
        @test node.id == 2
        @test node.time == 1.0
        @test node.parent == 1
        @test node.kind == SampledLeaf
        @test branch_length(tree, 2) == 1.0
        @test branch_length(tree, 3) == 1.2
        @test_throws ErrorException branch_length(tree, 1)
    end
end
