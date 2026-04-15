struct TreePlot
    svg::String
end

Base.show(io::IO, plot::TreePlot) = print(io, "TreePlot()")
Base.show(io::IO, ::MIME"image/svg+xml", plot::TreePlot) = print(io, plot.svg)

"""
    plot_tree(tree::Tree; width=800, height=500, node_color=nothing, labels=nothing,
              node_size=4.0, edge_width=1.4, edge_color=:default)

Return a lightweight SVG-backed `TreePlot` for quick tree inspection.

The plot is built from `node_positions(tree)` and `edge_segments(tree)`, using
the package's preorder tip/layout conventions. If `node_color` is supplied, it
must be node-index-aligned. Numeric colors are mapped continuously; other values
are mapped as discrete categories. If `labels` is supplied, it must also be
node-index-aligned.

This is a static inspection viewer, not a general plotting framework.
"""
function plot_tree(tree::Tree;
                   width::Real=800,
                   height::Real=500,
                   node_color=nothing,
                   labels=nothing,
                   node_size::Real=4.0,
                   edge_width::Real=1.4,
                   edge_color=:default)
    plot_width = _positive_float(width, "width")
    plot_height = _positive_float(height, "height")
    marker_size = _positive_float(node_size, "node_size")
    line_width = _positive_float(edge_width, "edge_width")
    line_color = _edge_color(edge_color)

    _check_aligned_vector(tree, node_color, "node_color"; reject_nothing=true)
    _check_aligned_vector(tree, labels, "labels"; reject_nothing=true)

    x, y = node_positions(tree)
    x1, y1, x2, y2 = edge_segments(tree)
    colors = _node_colors(tree, node_color)
    label_text = labels === nothing ? nothing : string.(labels)

    svg = _tree_svg(tree, x, y, x1, y1, x2, y2, colors, label_text, marker_size, line_width, line_color, plot_width, plot_height)
    return TreePlot(svg)
end

function _positive_float(value::Real, name::String)
    out = float(value)
    isfinite(out) && out > 0 || throw(ArgumentError("$name must be a positive finite value."))
    return out
end

function _check_aligned_vector(tree::Tree, values, name::String; reject_nothing::Bool=false)
    values === nothing && return nothing
    length(values) == length(tree) || throw(ArgumentError("$name must have length $(length(tree)); got $(length(values))."))

    if any(ismissing, values)
        throw(ArgumentError("$name must not contain missing values."))
    end
    if reject_nothing && any(value -> value === nothing, values)
        throw(ArgumentError("$name must not contain nothing values."))
    end

    return nothing
end

function _edge_color(edge_color)
    (edge_color === nothing || edge_color === :default) && return "#444444"
    edge_color === missing && throw(ArgumentError("edge_color must not be missing."))
    return string(edge_color)
end

function _node_colors(tree::Tree, node_color)
    node_color === nothing && return _default_node_colors(tree)
    isempty(node_color) && return String[]

    if all(value -> value isa Real, node_color)
        return _numeric_colors(Float64.(node_color))
    end

    return _categorical_colors(node_color)
end

function _default_node_colors(tree::Tree)
    colors = String[]
    root_set = Set(roots(tree))

    for i in eachindex(tree)
        if i in root_set
            push!(colors, "#d95f02")
        elseif isleaf(tree, i)
            push!(colors, "#1b9e77")
        else
            push!(colors, "#2f6f9f")
        end
    end

    return colors
end

function _numeric_colors(values::Vector{Float64})
    lo = minimum(values)
    hi = maximum(values)
    span = hi - lo

    if span == 0.0
        return fill("#2f6f9f", length(values))
    end

    return [_interpolate_color((value - lo) / span) for value in values]
end

function _interpolate_color(t::Float64)
    r1, g1, b1 = 42, 111, 159
    r2, g2, b2 = 217, 95, 2
    r = round(Int, r1 + t * (r2 - r1))
    g = round(Int, g1 + t * (g2 - g1))
    b = round(Int, b1 + t * (b2 - b1))
    return _hex_color(r, g, b)
end

function _categorical_colors(values)
    palette = ["#2f6f9f", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#666666"]
    seen = Dict{Any, String}()
    out = String[]

    for value in values
        color = get!(seen, value) do
            palette[mod1(length(seen) + 1, length(palette))]
        end
        push!(out, color)
    end

    return out
end

function _hex_color(r::Int, g::Int, b::Int)
    return "#" * string(r, base=16, pad=2) * string(g, base=16, pad=2) * string(b, base=16, pad=2)
end

function _tree_svg(tree::Tree, x, y, x1, y1, x2, y2, colors, labels, node_size::Float64, edge_width::Float64, edge_color::String, width::Float64, height::Float64)
    margin = 32.0

    if length(tree) == 0
        return """
        <svg xmlns="http://www.w3.org/2000/svg" width="$(round(Int, width))" height="$(round(Int, height))" viewBox="0 0 $(round(Int, width)) $(round(Int, height))">
          <rect width="100%" height="100%" fill="white"/>
        </svg>
        """
    end

    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    xspan = xmax - xmin
    yspan = ymax - ymin

    sx(value) = margin + (xspan == 0.0 ? 0.5 : (value - xmin) / xspan) * (width - 2margin)
    sy(value) = height - margin - (yspan == 0.0 ? 0.5 : (value - ymin) / yspan) * (height - 2margin)

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(round(Int, width))" height="$(round(Int, height))" viewBox="0 0 $(round(Int, width)) $(round(Int, height))">""")
    println(io, """  <rect width="100%" height="100%" fill="white"/>""")
    println(io, """  <g fill="none" stroke="$(_escape_xml(edge_color))" stroke-width="$edge_width" stroke-linecap="round">""")

    for k in eachindex(x1)
        println(io, """    <line x1="$(sx(x1[k]))" y1="$(sy(y1[k]))" x2="$(sx(x2[k]))" y2="$(sy(y2[k]))"/>""")
    end

    println(io, "  </g>")
    println(io, """  <g stroke="white" stroke-width="1">""")

    for i in eachindex(tree)
        println(io, """    <circle cx="$(sx(x[i]))" cy="$(sy(y[i]))" r="$node_size" fill="$(colors[i])"/>""")
    end

    println(io, "  </g>")

    if labels !== nothing
        println(io, """  <g font-family="sans-serif" font-size="11" fill="#222222">""")
        for i in eachindex(tree)
            println(io, """    <text x="$(sx(x[i]) + node_size + 3)" y="$(sy(y[i]) + 4)">$(_escape_xml(labels[i]))</text>""")
        end
        println(io, "  </g>")
    end

    println(io, "</svg>")
    return String(take!(io))
end

function _escape_xml(value::AbstractString)
    out = replace(value, "&" => "&amp;")
    out = replace(out, "<" => "&lt;")
    out = replace(out, ">" => "&gt;")
    out = replace(out, "\"" => "&quot;")
    return replace(out, "'" => "&apos;")
end
