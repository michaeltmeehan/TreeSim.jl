using Test
using TreeSim
using EpiSim

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

@testset "EpiSim EventLog adaptor" begin
    @testset "serial samples on sibling transmissions produce reconstructed tree" begin
        log = EventLog(
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [1, 2, 3, 2, 3],
            [0, 1, 1, 0, 0],
            [EK_Seeding, EK_Transmission, EK_Transmission, EK_SerialSampling, EK_SerialSampling],
        )

        tree = tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test nnodes(tree) == 4
        @test tree.time == [0.0, 1.0, 3.0, 4.0]
        @test tree.kind == [Root, Binary, SampledLeaf, SampledLeaf]
        @test tree.host == [1, 1, 2, 3]
        @test children(tree, 1) == [2]
        @test children(tree, 2) == [4, 3]
        @test all(!=(UnsampledUnary), tree.kind)
        @test sort(tree.label[tips(tree)]) == [2, 3]

        forest = forest_from_eventlog(log)
        @test length(forest) == 1
        @test forest[1].time == tree.time
        @test forest[1].host == tree.host
    end

    @testset "fossilized sampling becomes unary when sampled descendants remain" begin
        log = EventLog(
            [0.0, 1.0, 2.0, 3.0],
            [1, 2, 2, 2],
            [0, 1, 0, 0],
            [EK_Seeding, EK_Transmission, EK_FossilizedSampling, EK_SerialSampling],
        )

        tree = tree_from_eventlog(log)
        full = full_tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test tree.kind == [Root, SampledUnary, SampledLeaf]
        @test tree.host == [1, 2, 2]
        @test tree.label == [0, 0, 2]
        @test ancestors(tree, 3) == [2, 1]
        @test full.kind == [Root, UnsampledUnary, SampledUnary, SampledLeaf]
        @test validate_full_tree_against_eventlog(log, full)
        @test ancestors(full, 4) == [3, 2, 1]
    end

    @testset "unsampled side lineages are pruned" begin
        log = EventLog(
            [0.0, 1.0, 2.0, 3.0],
            [1, 2, 3, 2],
            [0, 1, 1, 0],
            [EK_Seeding, EK_Transmission, EK_Transmission, EK_SerialSampling],
        )

        tree = tree_from_eventlog(log)
        full = full_tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test tree.kind == [Root, SampledLeaf]
        @test tree.host == [1, 2]
        @test all(!=(3), tree.host)
        @test full.kind == [Root, UnsampledUnary, SampledLeaf]
        @test full.host == [1, 1, 2]
        @test validate_full_tree_against_eventlog(log, full)
    end

    @testset "reconstructed tree aliases collapse unsampled unary intermediates" begin
        log = EventLog(
            [0.0, 1.0, 2.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_SerialSampling],
        )

        tree = reconstructed_tree_from_eventlog(log)

        @test tree_from_eventlog(log).kind == tree.kind
        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test tree.time == [0.0, 2.0]
        @test tree.kind == [Root, SampledLeaf]
        @test tree.host == [1, 2]
        @test branch_length(tree, 2) == 2.0
        @test all(!=(UnsampledUnary), tree.kind)
    end

    @testset "full tree retains sampled-ancestry transmission intermediates" begin
        log = EventLog(
            [0.0, 1.0, 2.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_SerialSampling],
        )

        tree = full_tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_full_tree_against_eventlog(log, tree)
        @test tree.time == [0.0, 1.0, 2.0]
        @test tree.kind == [Root, UnsampledUnary, SampledLeaf]
        @test tree.host == [1, 1, 2]
        @test children(tree, 1) == [2]
        @test children(tree, 2) == [3]
        @test branch_length(tree, 2) == 1.0
        @test branch_length(tree, 3) == 1.0
    end

    @testset "reconstructed binary validation allows collapsed immediate infectee" begin
        log = EventLog(
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [1, 2, 1, 3, 3],
            [0, 1, 0, 2, 0],
            [EK_Seeding, EK_Transmission, EK_SerialSampling, EK_Transmission, EK_SerialSampling],
        )

        tree = tree_from_eventlog(log)
        full = full_tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test tree.kind == [Root, Binary, SampledLeaf, SampledLeaf]
        @test tree.host == [1, 1, 1, 3]
        @test all(!=(2), tree.host)
        @test_throws ErrorException validate_full_tree_against_eventlog(log, tree)

        @test validate_tree(full)
        @test validate_full_tree_against_eventlog(log, full)
        @test full.kind == [Root, Binary, SampledLeaf, UnsampledUnary, SampledLeaf]
        @test 2 in full.host
    end

    @testset "logs without samples produce an empty tree" begin
        log = EventLog(
            [0.0, 1.0, 2.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_Removal],
        )

        tree = tree_from_eventlog(log)
        full = full_tree_from_eventlog(log)

        @test isempty(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test validate_full_tree_against_eventlog(log, full)
        @test forest_from_eventlog(log) == Tree[]
    end

    @testset "bridge validation catches event-node mismatches" begin
        log = EventLog(
            [0.0, 1.0, 2.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_SerialSampling],
        )
        tree = tree_from_eventlog(log)
        shifted = Tree(
            copy(tree.time),
            copy(tree.left),
            copy(tree.right),
            copy(tree.parent),
            copy(tree.kind),
            [1, 1, 99],
            copy(tree.label),
        )

        @test_throws ErrorException validate_tree_against_eventlog(log, shifted)
    end

    @testset "multiple retained seeded components are rejected" begin
        log = EventLog(
            [0.0, 0.0, 1.0, 1.5],
            [1, 2, 1, 2],
            [0, 0, 0, 0],
            [EK_Seeding, EK_Seeding, EK_SerialSampling, EK_SerialSampling],
        )

        err = try
            tree_from_eventlog(log)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("strict single-tree API", sprint(showerror, err))

        forest = forest_from_eventlog(log)
        @test length(forest) == 2
        @test all(validate_tree, forest)
        @test validate_tree_against_eventlog(log, forest)
        @test [tree.host[root(tree)] for tree in forest] == [1, 2]
        @test [only(tree.label[tips(tree)]) for tree in forest] == [1, 2]
    end

    @testset "forest components are deterministically ordered by root time and host" begin
        log = EventLog(
            [0.0, 0.0, 1.0, 1.5],
            [5, 2, 5, 2],
            [0, 0, 0, 0],
            [EK_Seeding, EK_Seeding, EK_SerialSampling, EK_SerialSampling],
        )

        forest = forest_from_eventlog(log)

        @test length(forest) == 2
        @test [tree.time[root(tree)] for tree in forest] == [0.0, 0.0]
        @test [tree.host[root(tree)] for tree in forest] == [2, 5]
        @test [only(tree.label[tips(tree)]) for tree in forest] == [2, 5]
    end

    @testset "forest validation counts events globally" begin
        log = EventLog(
            [0.0, 0.0, 1.0, 1.5],
            [1, 2, 1, 2],
            [0, 0, 0, 0],
            [EK_Seeding, EK_Seeding, EK_SerialSampling, EK_SerialSampling],
        )
        forest = forest_from_eventlog(log)
        duplicated = [forest[1], forest[1]]

        err = try
            validate_tree_against_eventlog(log, duplicated)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("does not correspond to Seeding event", sprint(showerror, err))
    end

    @testset "equal-time retained full ancestry is rejected clearly" begin
        log = EventLog(
            [0.0, 0.0, 1.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_SerialSampling],
        )

        err = try
            full_tree_from_eventlog(log)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("Retained ancestry times must be strictly increasing", sprint(showerror, err))
    end

    @testset "duplicate retained sampling signatures are rejected" begin
        log = EventLog(
            [0.0, 1.0, 1.0],
            [1, 1, 1],
            [0, 0, 0],
            [EK_Seeding, EK_FossilizedSampling, EK_FossilizedSampling],
        )

        err = try
            tree_from_eventlog(log)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("Ambiguous duplicate retained sampling event signature", sprint(showerror, err))
    end

    @testset "sparse positive host ids are accepted" begin
        log = EventLog(
            [0.0, 2.0],
            [100, 100],
            [0, 0],
            [EK_Seeding, EK_SerialSampling],
        )

        tree = tree_from_eventlog(log)

        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test tree.host == [100, 100]
        @test tree.label == [0, 100]
    end
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

function nonzero_child_slots(tree)
    count = 0
    for i in eachindex(tree)
        count += tree.left[i] != 0
        count += tree.right[i] != 0
    end
    return count
end

function test_edge_segments_match_node_positions(tree)
    parents, child_nodes = parent_child_pairs(tree)
    x, y = node_positions(tree)
    x1, y1, x2, y2 = edge_segments(tree)

    @test length(parents) == length(child_nodes)
    @test length(x1) == length(parents)
    @test length(y1) == length(parents)
    @test length(x2) == length(parents)
    @test length(y2) == length(parents)

    for k in eachindex(parents)
        @test x1[k] == x[parents[k]]
        @test y1[k] == y[parents[k]]
        @test x2[k] == x[child_nodes[k]]
        @test y2[k] == y[child_nodes[k]]
    end
end

function svg_string(plot)
    return sprint(show, MIME("image/svg+xml"), plot)
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

        tip_y = tip_positions(tree)
        x, y = node_positions(tree)
        @test length(tip_y) == nnodes(tree)
        @test isnan(tip_y[1])
        @test isnan(tip_y[2])
        @test tip_y[tips(tree)] ≈ [1.0, 2.0, 3.0]
        @test x == node_depths(tree)
        @test y[tips(tree)] ≈ tip_y[tips(tree)]
        @test y[2] ≈ (y[4] + y[5]) / 2
        @test y[1] ≈ (y[2] + y[3]) / 2

        parents, child_nodes = parent_child_pairs(tree)
        @test parents == [1, 1, 2, 2]
        @test child_nodes == [2, 3, 4, 5]
        @test length(parents) == nonzero_child_slots(tree)
        x1, y1, x2, y2 = edge_segments(tree)
        @test x1 ≈ [0.0, 0.0, 0.4, 0.4]
        @test y1 ≈ [2.25, 2.25, 1.5, 1.5]
        @test x2 ≈ [0.4, 0.8, 1.1, 1.3]
        @test y2 ≈ [1.5, 3.0, 1.0, 2.0]
        test_edge_segments_match_node_positions(tree)
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

        tip_y = tip_positions(tree)
        x, y = node_positions(tree)
        @test all(isnan, tip_y[1:3])
        @test tip_y[4] == 1.0
        @test x == node_depths(tree)
        @test y[4] == 1.0
        @test y[1] == y[2] == y[3] == y[4]

        parents, child_nodes = parent_child_pairs(tree)
        @test parents == [1, 2, 3]
        @test child_nodes == [2, 3, 4]
        @test length(parents) == nonzero_child_slots(tree)
        x1, y1, x2, y2 = edge_segments(tree)
        @test x1 ≈ [0.0, 0.4, 0.9]
        @test y1 ≈ [1.0, 1.0, 1.0]
        @test x2 ≈ [0.4, 0.9, 1.5]
        @test y2 ≈ [1.0, 1.0, 1.0]
        test_edge_segments_match_node_positions(tree)
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

        tip_y = tip_positions(tree)
        x, y = node_positions(tree)
        @test tip_y[tips(tree)] ≈ [1.0, 2.0, 3.0, 4.0]
        @test all(isnan, tip_y[internal_nodes(tree)])
        @test x == node_depths(tree)
        @test y[2] ≈ 1.5
        @test y[3] ≈ 3.5
        @test y[1] ≈ 2.5

        parents, child_nodes = parent_child_pairs(tree)
        @test parents == [1, 1, 2, 2, 3, 3]
        @test child_nodes == [2, 3, 4, 5, 6, 7]
        @test length(parents) == nonzero_child_slots(tree)
        test_edge_segments_match_node_positions(tree)
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

    @testset "tree collection interface" begin
        tree = canonical_binary_tree()
        empty = Tree()

        @test !isempty(tree)
        @test isempty(empty)
        @test firstindex(tree) == 1
        @test lastindex(tree) == nnodes(tree)
        @test firstindex(empty) == 1
        @test lastindex(empty) == 0
        @test eltype(Tree) == Node
        @test Base.IteratorSize(Tree) == Base.HasLength()
        @test Base.IteratorEltype(Tree) == Base.HasEltype()

        @test tree[end] == tree[nnodes(tree)]
        @test first(tree) == tree[1]
        @test last(tree) == tree[end]

        nodes = collect(tree)
        @test nodes isa Vector{Node}
        @test [node.id for node in nodes] == collect(eachindex(tree))
        @test nodes == [tree[i] for i in eachindex(tree)]

        @test_throws BoundsError tree[0]
        @test_throws BoundsError tree[nnodes(tree) + 1]
        @test_throws BoundsError empty[end]
    end

    @testset "minimal plotting" begin
        binary = canonical_binary_tree()
        unary = canonical_unary_tree()

        empty_plot = plot_tree(Tree())
        binary_plot = plot_tree(binary)
        unary_plot = plot_tree(unary; node_color=collect(1:nnodes(unary)))
        sized_plot = plot_tree(binary; width=640, height=360, edge_color="black")

        @test occursin("<svg", svg_string(empty_plot))
        @test occursin("<svg", svg_string(binary_plot))
        @test occursin("<line", svg_string(binary_plot))
        @test occursin("<circle", svg_string(binary_plot))
        @test occursin("#d95f02", svg_string(binary_plot))
        @test occursin("#1b9e77", svg_string(binary_plot))
        @test occursin("<svg", svg_string(unary_plot))
        @test occursin("width=\"640\"", svg_string(sized_plot))
        @test occursin("height=\"360\"", svg_string(sized_plot))
        @test occursin("stroke=\"black\"", svg_string(sized_plot))

        categorical_plot = plot_tree(binary; node_color=binary.kind)
        @test occursin("<circle", svg_string(categorical_plot))

        labeled_plot = plot_tree(binary; labels=collect(eachindex(binary)))
        string_labeled_plot = plot_tree(binary; labels=["n$i" for i in eachindex(binary)])
        @test occursin(">1<", svg_string(labeled_plot))
        @test occursin("n1", svg_string(string_labeled_plot))

        @test_throws ArgumentError plot_tree(binary; width=0)
        @test_throws ArgumentError plot_tree(binary; height=-1)
        @test_throws ArgumentError plot_tree(binary; node_color=[1, 2])
        @test_throws ArgumentError plot_tree(binary; node_color=[1, 2, missing, 4, 5])
        @test_throws ArgumentError plot_tree(binary; node_color=Any[1, 2, nothing, 4, 5])
        @test_throws ArgumentError plot_tree(binary; labels=["too", "short"])
        @test_throws ArgumentError plot_tree(binary; labels=Any["n1", "n2", missing, "n4", "n5"])
        @test_throws ArgumentError plot_tree(binary; labels=Any["n1", "n2", nothing, "n4", "n5"])
        @test_throws ArgumentError plot_tree(binary; edge_color=missing)
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
        @test tip_positions(empty) == Float64[]
        x_empty, y_empty = node_positions(empty)
        @test x_empty == Float64[]
        @test y_empty == Float64[]
        parents_empty, children_empty = parent_child_pairs(empty)
        @test parents_empty == Int[]
        @test children_empty == Int[]
        x1_empty, y1_empty, x2_empty, y2_empty = edge_segments(empty)
        @test x1_empty == Float64[]
        @test y1_empty == Float64[]
        @test x2_empty == Float64[]
        @test y2_empty == Float64[]

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

        tip_y = tip_positions(multi_root_tree)
        x, y = node_positions(multi_root_tree)
        @test all(isnan, tip_y[1:2])
        @test tip_y[3:4] == [1.0, 2.0]
        @test tip_y[tips(multi_root_tree)] ≈ [1.0, 2.0]
        @test x == node_depths(multi_root_tree)
        @test y == [1.0, 2.0, 1.0, 2.0]

        parents, child_nodes = parent_child_pairs(multi_root_tree)
        @test parents == [1, 2]
        @test child_nodes == [3, 4]
        @test length(parents) == nonzero_child_slots(multi_root_tree)
        x1, y1, x2, y2 = edge_segments(multi_root_tree)
        @test x1 ≈ [0.0, 0.0]
        @test y1 ≈ [1.0, 2.0]
        @test x2 ≈ [0.5, 0.7]
        @test y2 ≈ [1.0, 2.0]
        test_edge_segments_match_node_positions(multi_root_tree)

        multi_root_plot = plot_tree(multi_root_tree; node_color=multi_root_tree.kind)
        @test occursin("<svg", svg_string(multi_root_plot))
    end

    @testset "node row view and branch length" begin
        tree = canonical_binary_tree()
        root_node = tree[1]
        internal_node = tree[2]
        tip_node = tree[4]

        @test tip_node isa Node
        @test tip_node.id == 4
        @test tip_node.time == 1.1
        @test tip_node.parent == 2
        @test tip_node.kind == SampledLeaf
        @test branch_length(tree, 2) == 0.4
        @test branch_length(tree, 4) ≈ 0.7
        @test_throws ErrorException branch_length(tree, 1)

        tree_display = sprint(show, tree)
        @test startswith(tree_display, "Tree(")
        @test occursin("5 nodes", tree_display)
        @test occursin("1 root", tree_display)
        @test occursin("2 internal", tree_display)
        @test occursin("3 tips", tree_display)

        root_display = sprint(show, root_node)
        @test occursin("Node(1, Root", root_display)
        @test occursin("time=0.0", root_display)
        @test occursin("children=[2, 3]", root_display)
        @test !occursin("parent=0", root_display)

        internal_display = sprint(show, internal_node)
        @test occursin("Node(2, Binary", internal_display)
        @test occursin("parent=1", internal_display)
        @test occursin("children=[4, 5]", internal_display)

        tip_display = sprint(show, tip_node)
        @test occursin("Node(4, SampledLeaf", tip_display)
        @test occursin("time=1.1", tip_display)
        @test occursin("parent=2", tip_display)
        @test occursin("label=102", tip_display)
        @test !occursin("left=0", tip_display)
        @test !occursin("right=0", tip_display)
        @test !occursin("children=", tip_display)
        @test !occursin("host=0", tip_display)

        metadata_node = Node(10, 2.5, 0, 0, 4, SampledLeaf, 7, 9001)
        metadata_display = sprint(show, metadata_node)
        @test occursin("host=7", metadata_display)
        @test occursin("label=9001", metadata_display)
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
