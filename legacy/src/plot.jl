@recipe function f(phylo::Phylogeny)

    timespan = extrema([node.data["t"] for node in traversal(phylo.tree, inorder)])
    xlims := timespan
    xscale := :identity
    legend := nothing
    xlabel := "Time"
    yaxis := false

    @series begin
        linewidth := 5
        # linecolor := [node.data["host"] for node in traversal(phylo.tree, preorder)]
        # TODO: Map host IDs to colors
        linecolor := reverse([branch.data["host"] for branch in getbranches(phylo.tree)])

        markersize := 10
        # TODO: Map host IDs to colors
        markercolor := [node.id for node in traversal(phylo.tree, preorder)]
        markerstrokewidth := 0
        seriestype := :path
        phylo.tree
    end

end