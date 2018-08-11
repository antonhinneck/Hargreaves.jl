__precompile__(true)
module Hargreaves
export wireplot, gridplot

include("layout.jl")
using Cairo, Colors, LightGraphs

function _get_weight(weight_matrix::Array{T,2} where T <: Number, e::AbstractEdge{T} where T <: Integer)
#This function only supports upper triangular matrices.
#(only unidirected graphs)

    weigth_matrix = UpperTriangular(weight_matrix)
    weight = 0
    if src(e) < dst(e) && src(e) != dst(e)
        weight = weight_matrix[src(e), dst(e)]
    elseif dst(e) < src(e) && src(e) != dst(e)
        weight = weight_matrix[dst(e), src(e)]
    end
    return weight
end

function _get_width(value::Number,
                    actual_min::Number,
                    actual_max::Number,
                    min_width::Number,
                    max_width::Number,
                    static_width_value::Number,
                    width_scale::Symbol,
                    static_widths::Bool)

    output_width = static_width_value
    normrange = max_width - min_width
    weight_range = actual_max - actual_min

    if value == 0
        output_width = 0
    else
        if !static_widths && !(normrange == 0) && !(weight_range == 0) && width_scale == :linear
            output_width = (value - actual_min) / (weight_range) * normrange
        elseif !static_widths && !(normrange == 0) && !(weight_range == 0) && width_scale == :quadratic
            # A proper quadratic function would be better.
            # The value max_width is exceeded by min_width * (actual_min / value) ^ 8.
            output_width = ((value - actual_min) / weight_range) ^ 2 * normrange + min_width * (actual_min / value) ^ 8
        end
        return output_width
    end
end

function _compute_node_positions(res_x, res_y, node_count, ratio=0.8, switch=1)
    center = [res_x / 2, res_y / 2]
    degrees = 2 * pi
    degree_sections = degrees / node_count
    main_radius = res_x / 2 * ratio
    positions = Array{Float64}(node_count, 2)
    for i in 1:node_count
        positions[i,1] = cos(degree_sections * (i - 1)) * main_radius + center[1]
        positions[i,2] = sin(degree_sections * (i - 1)) * main_radius + center[2]
    end
    return positions
end

function _get_color(value::Number,
                    actual_min::Number,
                    actual_max::Number,
                    static_colors::Bool,
                    color_scale,
                    scheme::Symbol = :generic)

    # This function only supports decreasing color values
    # while weights increase: IF-clause needed!

    (r_min, r_max) = (200, 112)
    (g_min, g_max) = (255, 30)
    (b_min, b_max) = (225, 112)
    if scheme == :sky
        (r_min, r_max) = (254, 8)
        (g_min, g_max) = (255, 29)
        (b_min, b_max) = (217, 88)
    elseif scheme == :blue
        (r_min, r_max) = (230, 0)
        (g_min, g_max) = (238, 55)
        (b_min, b_max) = (255, 160)
    elseif scheme == :purple
        (r_min, r_max) = (255, 75)
        (g_min, g_max) = (246, 0)
        (b_min, b_max) = (242, 106)
    end

    weight_range = actual_max - actual_min

    @inline function output_color(color_min, color_max)

        output = color_min / 255
        color_range = abs(color_max - color_min)

        if color_scale == :linear
            output = color_min
            if color_max < color_min
                output = (color_min - ((color_range) / weight_range) * value) / 255
            elseif color_max > color_min
                output = (color_max - ((color_range) / weight_range) * value) / 255
            end
        elseif color_scale == :quadratic
            output = color_min
            if color_max < color_min
                output = ((color_range / weight_range) * value) / 255
                output = color_min / 255 - output ^ 2
            elseif color_max > color_min
                output = (((color_range / weight_range) * value) ^ 2) / 255
            end
        end

        return output
    end # inline function

    r_out = output_color(r_min, r_max)
    g_out = output_color(g_min, g_max)
    b_out = output_color(b_min, b_max)

    if static_colors
    # Override min color used in legend's gradient
        r_min = r_max
        g_min = g_max
        b_min = b_max
    end
    return [
        [r_out,g_out,b_out],
        [0.0,0.0,0.0],
        [r_min / 255,g_min / 255,b_min / 255],
        [r_max / 255,g_max / 255,b_max / 255]
    ]
end

function wireplot(g::AbstractGraph{T=Int64}, basefn = "wireplot";
    distmx = weights(g),
    res_x = 4096, res_y = 4096,
    plot_to_surface_ratio = 0.9,
    wireplot_bg_color = [1.0,1.0,1.0],
    wireplot_jitter = 0,
    wireplot_node_diameter = 8.0,
    wireplot_node_label_offset = 16.0,
    wireplot_node_label_font_color = [0.0,0.0,0.0],
    wireplot_node_label_font_size = 48.0,
    wireplot_edge_color_scheme =:generic,
    wireplot_legend_font_color = [0.0,0.0,0.0],
    wireplot_legend_heading = "Edges' weights",
    wireplot_legend_show_values = true,
    static_widths = false,
    static_colors = false,
    color_scale = :linear,
    width_scale = :linear,
    min_width = 1,  # minimum normalized arc width
    max_width = 8, # maximum normalized arc width
    static_width_value = 1, # static arc width (if static_widths)
    export_type = :svg
    )

    (actual_min, actual_max) = extrema(distmx) #Get minimal and maximal weights

    @inline function _drop_zeros(array::Array{T1, 2} where T1 <: Number, max::T2 where T2 <: Number)
        # Drop weights if weight <= 0.
        output_min = max
        for a in array
            if a > 0 && a < output_min

                output_min = a

            end
        end
        return output_min
    end

    actual_min = _drop_zeros(distmx, actual_max)
    #println(actual_min)

    ## SETUP SURFACE
    #-----------------------

    center = [res_x / 2,res_y / 2]

    if export_type == :pdf
        c = CairoPDFSurface("$(basefn).pdf", res_x, res_y)
    else
        c = CairoSVGSurface("$(basefn).svg", res_x, res_y)
    end
    cr = CairoContext(c)
    select_font_face(cr, "Times", 1, 1)
    set_font_size(cr, wireplot_node_label_font_size)

    ## DRAW BACKGROUND
    #-----------------------

    set_source_rgb(cr, wireplot_bg_color...)
    rectangle(cr, 0.0, 0.0, res_x, res_y)
    fill(cr)

    ## PLOT EDGES
    #-----------------------

    plot_border_top = res_y
    plot_border_bottom = 0
    plot_border_left = res_x
    plot_border_right = 0

    node_pos = _compute_node_positions(res_x, res_y, nv(g), plot_to_surface_ratio, 1)

    degs = degree(g)

    for (i, e) in enumerate(edges(g))
        (s, d) = (src(e), dst(e))
        # indeg_d = indegs[d]
        # deg_d = indeg_d + outdeg_d

        if distmx == weights(g) || static_colors
            set_source_rgb(cr, _get_color(actual_max, actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[1]...)
        else
            set_source_rgb(cr, _get_color(_get_weight(distmx, e), actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[1]...)
        end
        line_width = _get_width(_get_weight(distmx, e),
                                      actual_min,
                                      actual_max,
                                      min_width,
                                      max_width,
                                      static_width_value,
                                      width_scale,
                                      static_widths)
        set_line_width(cr, line_width)
        if line_width > 0
            curve_to(cr, node_pos[s, 1], node_pos[s, 2],
                center[1] + rand() * wireplot_jitter - rand() * wireplot_jitter,
                center[2] + rand() * wireplot_jitter - rand() * wireplot_jitter,
                node_pos[d, 1], node_pos[d, 2])
            stroke(cr)
        end
    end
    set_line_width(cr, 1)

    ## PLOT NODES BASED ON COMPUTED POSITIONS
    #----------------------------------------

    check = false

    for i in vertices(g)

        (pos_x, pos_y) = node_pos[i,:]
        set_source_rgb(cr, [0.0,0.0,0.0]...)
        move_to(cr, pos_x, pos_y) # Prevents artifacts in the exported pdf file
        circle(cr, pos_x, pos_y, wireplot_node_diameter)
        fill(cr)
        text = string(i)
        label_extents = text_extents(cr, text)

        set_source_rgb(cr, wireplot_node_label_font_color...)

        if round(Int, pos_x) == center[1] && pos_y > center[2]
            label_origin_x = pos_x - label_extents[3] / 2 - 1
            label_origin_y = pos_y + wireplot_node_diameter / 2 + label_extents[4] + wireplot_node_label_offset
            label_border_left = pos_x - label_extents[3] / 2 - 1
            label_border_right = label_origin_x + label_extents[3]
            label_border_top = label_origin_y + label_extents[4]
            label_border_bottom = label_origin_y
        elseif round(Int, pos_x) == center[1] && pos_y < center[2]
            label_origin_x = pos_x - label_extents[3] / 2 - 1
            label_origin_y = pos_y - wireplot_node_diameter / 2 - wireplot_node_label_offset
            label_border_left = pos_x - label_extents[3] / 2 - 1
            label_border_right = label_origin_x + label_extents[3]
            label_border_top = label_origin_y + label_extents[4]
            label_border_bottom = label_origin_y
        elseif pos_x < center[1]
            label_origin_x = pos_x - wireplot_node_diameter / 2 - label_extents[3] - wireplot_node_label_offset
            label_origin_y = pos_y + label_extents[4] / 2 + sin(((i - 1) * 2 * pi) / nv(g)) * (wireplot_node_label_offset + wireplot_node_diameter / 2 + label_extents[4] / 2)
            label_border_left = label_origin_x
            label_border_right = label_origin_x + label_extents[3]
            label_border_top = label_origin_y + label_extents[4]
            label_border_bottom = label_origin_y
        elseif pos_x > center[1]
            label_origin_x = pos_x + wireplot_node_diameter / 2 + wireplot_node_diameter / 2 + wireplot_node_label_offset
            label_origin_y = pos_y + label_extents[4] / 2 + sin(((i - 1) * 2 * pi) / nv(g)) * (wireplot_node_label_offset + wireplot_node_diameter / 2 + label_extents[4] / 2)
            label_border_left = label_origin_x
            label_border_right = label_origin_x + label_extents[3]
            label_border_top = label_origin_y + label_extents[4]
            label_border_bottom = label_origin_y
        end

        move_to(cr, label_origin_x, label_origin_y)
        show_text(cr, text)
        plot_border_bottom = max(plot_border_bottom, label_border_bottom)
        plot_border_right = max(plot_border_right, label_border_right)
        plot_border_top = min(plot_border_top, label_border_top)
        plot_border_left = min(plot_border_left, label_border_left)
    end
    ## PLOT LEGEND
    #-----------------------

    ## HEADING

    heading_extents = text_extents(cr, wireplot_legend_heading)

    set_source_rgb(cr, wireplot_legend_font_color...)
    move_to(cr, plot_border_right - res_x * 0.25, plot_border_bottom)
    show_text(cr, wireplot_legend_heading)
    line_to(cr, plot_border_right, plot_border_bottom)

    stroke(cr)

    ## COLOR DESCRIPTOR

    rectangle(cr, plot_border_right - res_x * 0.23 + heading_extents[3],  plot_border_bottom - res_y * 0.015,  res_x * 0.23 - heading_extents[3], res_y * 0.01)

    gradient = pattern_create_linear(plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom,  plot_border_right,  plot_border_bottom)

    pattern_add_color_stop_rgb(gradient, 0, _get_color(0, actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, _get_color(0, actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[4]...)
    set_source(cr, gradient)

    fill_preserve(cr)
    destroy(gradient)

    ## COLOR DESCRIPTOR FONT

    if wireplot_legend_show_values
        set_source_rgb(cr, wireplot_legend_font_color...)
        text_min = string("[", 1, ", ...")
        text_max = string(round(Int, actual_max), "]")

        move_to(cr, plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom - res_y * 0.02)
        show_text(cr, text_min)

        move_to(cr, plot_border_right - text_extents(cr, text_max)[3] - text_extents(cr, text_max)[1], plot_border_bottom - res_y * 0.02)
        show_text(cr, text_max)
    end

    ## EXPORT SVG AND PNG
    #--------------------

    if export_type == :png
        write_to_png(c, "$(basefn).png")
    else
        finish(c)
    end
end

function gridplot(g::AbstractGraph{T=Int64}, basefn = "wireplot";
    distmx = weights(g),
    res_x = 4096, res_y = 4096,
    plot_to_surface_ratio = 0.9,
    wireplot_bg_color = [1.0,1.0,1.0],
    wireplot_jitter = 0,
    wireplot_node_diameter = 8.0,
    wireplot_node_label_offset = 16.0,
    wireplot_node_label_font_color = [0.0,0.0,0.0],
    wireplot_node_label_font_size = 48.0,
    wireplot_edge_color_scheme =:generic,
    wireplot_legend_font_color = [0.0,0.0,0.0],
    wireplot_legend_heading = "Edges' weights",
    wireplot_legend_show_values = true,
    static_widths = false,
    static_colors = false,
    color_scale = :linear,
    width_scale = :linear,
    min_width = 1,  # minimum normalized arc width
    max_width = 8, # maximum normalized arc width
    static_width_value = 1, # static arc width (if static_widths)
    export_type = :svg)

    (actual_min, actual_max) = extrema(distmx) #Get minimal and maximal weights

    @inline function _drop_zeros(array::Array{T1, 2} where T1 <: Number, max::T2 where T2 <: Number)
        # Drop weights if weight <= 0.
        output_min = max
        for a in array
            if a > 0 && a < output_min

                output_min = a

            end
        end
        return output_min
    end

    actual_min = _drop_zeros(distmx, actual_max)
    #println(actual_min)

    ## SETUP SURFACE
    #-----------------------

    plot_layout = layout([convert(Int64, res_x), convert(Int64, res_y)])

    if export_type == :pdf
        c = CairoPDFSurface("$(basefn).pdf", res_x, res_y)
    else
        c = CairoSVGSurface("$(basefn).svg", res_x, res_y)
    end
    cr = CairoContext(c)
    select_font_face(cr, "Times", 1, 1)
    set_font_size(cr, wireplot_node_label_font_size)

    ## DRAW BACKGROUND
    #-----------------------

    set_source_rgb(cr, wireplot_bg_color...)
    rectangle(cr, 0.0, 0.0, plot_layout.canvas_dimensions[1], plot_layout.canvas_dimensions[2])
    fill(cr)

    ## PLOT GRID
    #-----------

    current_width = plot_layout.plot_pane_def[3] / length(distmx[:, 1])
    current_height = plot_layout.plot_pane_def[4] / length(distmx[1, :])

    for i in 1:length(distmx[:, 1])
        for j in 1:length(distmx[1, :])

            current_x = plot_layout.plot_pane_def[1] + plot_layout.plot_pane_def[3] / length(distmx[:, 1]) * (i - 1)
            current_y = plot_layout.plot_pane_def[2] + plot_layout.plot_pane_def[4] / length(distmx[1, :]) * (j - 1)

            rectangle(cr, current_x, current_y, current_width, current_height)
            set_source_rgb(cr, _get_color(distmx[i, j], actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[1]...)
            fill(cr)

        end
    end

    ## PLOT LEGEND
    #-----------------------

    ## HEADING

    heading_extents = text_extents(cr, wireplot_legend_heading)

    set_source_rgb(cr, wireplot_legend_font_color...)
    move_to(cr, plot_layout.legend_def[1], plot_layout.legend_def[2] + heading_extents[4])
    show_text(cr, wireplot_legend_heading)
    line_to(cr, plot_layout.legend_def[1] + plot_layout.legend_def[3], plot_layout.legend_def[2] + heading_extents[4])

    stroke(cr)
    #=
    ## COLOR DESCRIPTOR

    rectangle(cr, plot_border_right - res_x * 0.23 + heading_extents[3],  plot_border_bottom - res_y * 0.015,  res_x * 0.23 - heading_extents[3], res_y * 0.01)

    gradient = pattern_create_linear(plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom,  plot_border_right,  plot_border_bottom)

    pattern_add_color_stop_rgb(gradient, 0, _get_color(0, actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, _get_color(0, actual_min, actual_max, static_colors, color_scale, wireplot_edge_color_scheme)[4]...)
    set_source(cr, gradient)

    fill_preserve(cr)
    destroy(gradient)

    ## COLOR DESCRIPTOR FONT

    if wireplot_legend_show_values
        set_source_rgb(cr, wireplot_legend_font_color...)
        text_min = string("[", 1, ", ...")
        text_max = string(round(Int, actual_max), "]")

        move_to(cr, plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom - res_y * 0.02)
        show_text(cr, text_min)

        move_to(cr, plot_border_right - text_extents(cr, text_max)[3] - text_extents(cr, text_max)[1], plot_border_bottom - res_y * 0.02)
        show_text(cr, text_max)
    end
    =#
    ## EXPORT SVG AND PNG
    #--------------------

    if export_type == :png
        write_to_png(c, "$(basefn).png")
    else
        finish(c)
    end
end

end # module
