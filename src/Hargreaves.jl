module Hargreaves

using Cairo, Colors, DataFrames, LightGraphs

include("style.jl")
include("functions.jl")
include("arcs.jl")

function wireplot(g::AbstractGraph, res_x=4096, res_y=4096, plot_to_surface_ratio=0.9, plot_name="wireplot.svg")
    max_cons = Δ(g)
    center = [res_x / 2,res_y / 2]

    c = CairoSVGSurface("wireplot.svg", res_x, res_y)
    cr = CairoContext(c)

    select_font_face(cr, "Times", 1, 1)
    set_font_size(cr, 56.0)
    
    ## DRAW BACKGROUND
    #-----------------

    set_source_rgb(cr, wireplot_bg_color...)
    rectangle(cr, 0.0, 0.0, res_x, res_y)
    fill(cr)

    ## PLOT NODE CONNECTIONS
    #-----------------------

    plot_border_top = res_y
    plot_border_bottom = 0
    plot_border_left = res_x
    plot_border_right = 0
    
    node_pos = compute_node_positions(res_x, res_y, nv(g), plot_to_surface_ratio, 1)
    degs = degree(g)

    for (i, e) in enumerate(edges(g))
        (s, d) = (src(e), dst(e))
        # indeg_d = indegs[d]
        # deg_d = indeg_d + outdeg_d
        set_source_rgb(cr, get_color(degs[d], max_cons)[1]...)
        curve_to(cr, node_pos[s, 1], node_pos[s, 2],
            center[1] + rand() * wireplot_jitter - rand() * wireplot_jitter,
            center[2] + rand() * wireplot_jitter - rand() * wireplot_jitter,
            node_pos[d, 1], node_pos[d, 2])
        stroke(cr)
    end

    ## PLOT NODES BASED ON COMPUTED POSITIONS
    #----------------------------------------

    for i in vertices(g)
        (pos_x, pos_y) = node_pos[i,:]
        
        set_source_rgb(cr, get_color(degs[i], max_cons)[2]...)
        # line 200
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
            label_border_right = label_origin_x[i] + label_extents[3]
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
    #-------------

    ## HEADING

    heading = "Arc frequency"
    heading_extents = text_extents(cr, heading)

    set_source_rgb(cr, wireplot_legend_font_color...)
    move_to(cr, plot_border_right - res_x * 0.25, plot_border_bottom)
    show_text(cr, heading)
    line_to(cr, plot_border_right, plot_border_bottom)

    stroke(cr)

    ## COLOR DESCRIPTOR

    rectangle(cr, plot_border_right - res_x * 0.23 + heading_extents[3],  plot_border_bottom - res_y * 0.015,  res_x * 0.23 - heading_extents[3], res_y * 0.01)

    gradient = pattern_create_linear(plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom,  plot_border_right,  plot_border_bottom)

    pattern_add_color_stop_rgb(gradient, 0, get_color(0, max_cons)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, get_color(max_cons, max_cons)[4]...)
    set_source(cr, gradient)
    fill_preserve(cr)
    destroy(gradient)

    ## COLOR DESCRIPTOR FONT
    set_source_rgb(cr, wireplot_legend_font_color...)
    text_min = string("[", 1, ", ...")
    text_max = string(max_cons, "]")

    move_to(cr, plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom - res_y * 0.02)
    show_text(cr, text_min)

    move_to(cr, plot_border_right - text_extents(cr, text_max)[3] - text_extents(cr, text_max)[1], plot_border_bottom - res_y * 0.02)
    show_text(cr, text_max)

    ## EXPORT SVG AND PNG
    #--------------------

    write_to_png(c, "wireplot.png")

    finish(c)
        
end

#################
#################
#################

function wireplot(connections::Array{Int64,2}, res_x=4096, res_y=4096, plot_to_surface_ratio=0.9, plot_name="wireplot.svg")

    data = arcs(connections)

    center = [res_x / 2,res_y / 2]

    c = CairoSVGSurface("wireplot.svg", res_x, res_y)
    cr = CairoContext(c)

    select_font_face(cr, "Times", 1, 1)
    set_font_size(cr, 56.0)

    ## DRAW BACKGROUND
    #-----------------

    set_source_rgb(cr, wireplot_bg_color...)
    rectangle(cr, 0.0, 0.0, res_x, res_y)
    fill(cr)

    ## COMPUTE NODE POSITIONS
    #------------------------

    node_positions = compute_node_positions(res_x, res_y, length(data.arc_repo[:,1]), plot_to_surface_ratio, 1)

    ## PLOT NODE CONNECTIONS
    #-----------------------

    for i in 1:length(data.arcs_df[:,1])

        if data.arcs_df[i,1] > 0 || data.arcs_df[i,2] > 0

            set_source_rgb(cr, get_color(data.arcs_df[i,3], data.max_cons)[1]...)

            curve_to(cr, node_positions[data.arcs_df[i,1], 1],
                         node_positions[data.arcs_df[i,1], 2],
                         center[1] + rand() * wireplot_jitter - rand() * wireplot_jitter,
                         center[2] + rand() * wireplot_jitter - rand() * wireplot_jitter,
                         node_positions[data.arcs_df[i,2], 1],
                         node_positions[data.arcs_df[i,2], 2])
            stroke(cr)    
        end

    end

    ## PLOT NODES BASED ON COMPUTED POSITIONS
    #----------------------------------------

    for i in 1:length(data.arc_repo[:,1])

        set_source_rgb(cr, get_color(data.arcs_df[1,3], data.max_cons)[2]...)
        arc(cr, node_positions[i,1], node_positions[i,2], wireplot_node_diameter, wireplot_node_diameter, 2 * wireplot_node_diameter * pi)
        fill(cr)

    end

    ## PLOT NODES' LABELS
    #--------------------

    label_origin_x = Array{Float64,1}(length(node_positions[:,1]))
    label_origin_y = Array{Float64,1}(length(node_positions[:,1]))
    label_border_left = Array{Float64,1}(length(node_positions[:,1]))
    label_border_right = Array{Float64,1}(length(node_positions[:,1]))
    label_border_top = Array{Float64,1}(length(node_positions[:,1]))
    label_border_bottom = Array{Float64,1}(length(node_positions[:,1]))

    for i in 1:length(data.arc_repo[:,1])

        text = string(i)
        label_extents = text_extents(cr, text)

        set_source_rgb(cr, wireplot_node_label_font_color...)

        if round(Int, node_positions[i,1]) == center[1] && node_positions[i,2] > center[2]

            label_origin_x[i] = node_positions[i,1] - label_extents[3] / 2 - 1
            label_origin_y[i] = node_positions[i,2] + wireplot_node_diameter / 2 + label_extents[4] + wireplot_node_label_offset
            label_border_left[i] = node_positions[i,1] - label_extents[3] / 2 - 1
            label_border_right[i] = label_origin_x[i] + label_extents[3]
            label_border_top[i] = label_origin_y[i] + label_extents[4]
            label_border_bottom[i] = label_origin_y[i]

            move_to(cr, label_origin_x[i], label_origin_y[i])
            show_text(cr, string(i))

        elseif round(Int, node_positions[i,1]) == center[1] && node_positions[i,2] < center[2]

            label_origin_x[i] = node_positions[i,1] - label_extents[3] / 2 - 1
            label_origin_y[i] = node_positions[i,2] - wireplot_node_diameter / 2 - wireplot_node_label_offset
            label_border_left[i] = node_positions[i,1] - label_extents[3] / 2 - 1
            label_border_right[i] = label_origin_x[i] + label_extents[3]
            label_border_top[i] = label_origin_y[i] + label_extents[4]
            label_border_bottom[i] = label_origin_y[i]

            move_to(cr, label_origin_x[i], label_origin_y[i])
            show_text(cr, string(i))

        elseif node_positions[i,1] < center[1]

            label_origin_x[i] = node_positions[i,1] - wireplot_node_diameter / 2 - label_extents[3] - wireplot_node_label_offset
            label_origin_y[i] = node_positions[i,2] + label_extents[4] / 2 + sin(((i - 1) * 2 * pi) / length(node_positions[:,1])) * (wireplot_node_label_offset + wireplot_node_diameter / 2 + label_extents[4] / 2)
            label_border_left[i] = label_origin_x[i]
            label_border_right[i] = label_origin_x[i] + label_extents[3]
            label_border_top[i] = label_origin_y[i] + label_extents[4]
            label_border_bottom[i] = label_origin_y[i]

            move_to(cr, label_origin_x[i], label_origin_y[i])
            show_text(cr, string(i))

        elseif node_positions[i,1] > center[1]

            label_origin_x[i] = node_positions[i,1] + wireplot_node_diameter / 2 + wireplot_node_diameter / 2 + wireplot_node_label_offset
            label_origin_y[i] = node_positions[i,2] + label_extents[4] / 2 + sin(((i - 1) * 2 * pi) / length(node_positions[:,1])) * (wireplot_node_label_offset + wireplot_node_diameter / 2 + label_extents[4] / 2)
            label_border_left[i] = label_origin_x[i]
            label_border_right[i] = label_origin_x[i] + label_extents[3]
            label_border_top[i] = label_origin_y[i] + label_extents[4]
            label_border_bottom[i] = label_origin_y[i]

            move_to(cr, label_origin_x[i], label_origin_y[i])
            show_text(cr, string(i))

        end
    end

    fill(cr)

    ## COMPUTE PLOT BORDERS
    #----------------------

    plot_border_right = 0.0
    plot_border_top = res_y
    plot_border_left = res_x
    plot_border_bottom = 0.0

    for i in 1:length(node_positions[:,1])

        if label_border_right[i] > plot_border_right

            plot_border_right = label_border_right[i]

        end

        if label_border_top[i] < plot_border_top

            plot_border_top = label_border_top[i]

        end

        if label_border_left[i] < plot_border_left

            plot_border_left = label_border_left[i]

        end

        if label_border_bottom[i] > plot_border_bottom

            plot_border_bottom = label_border_bottom[i]

        end
    end

    ## PLOT LEGEND
    #-------------

        ## HEADING

    heading = "Arc frequency"
    heading_extents = text_extents(cr, heading)

    set_source_rgb(cr, wireplot_legend_font_color...)
    move_to(cr, plot_border_right - res_x * 0.25, plot_border_bottom)
    show_text(cr, heading)
    line_to(cr, plot_border_right, plot_border_bottom)

    stroke(cr)

        ## COLOR DESCRIPTOR

    rectangle(cr, plot_border_right - res_x * 0.23 + heading_extents[3],  plot_border_bottom - res_y * 0.015,  res_x * 0.23 - heading_extents[3], res_y * 0.01)

    gradient = pattern_create_linear(plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom,  plot_border_right,  plot_border_bottom)

    pattern_add_color_stop_rgb(gradient, 0, get_color(0, data.max_cons)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, get_color(data.max_cons, data.max_cons)[4]...)

    set_source(cr, gradient)

    fill_preserve(cr)

    destroy(gradient)

        ## COLOR DESCRIPTOR FONT

    set_source_rgb(cr, wireplot_legend_font_color...)

    text_min = string("[", 1, ", ...")
    text_max = string(data.max_cons, "]")

    move_to(cr, plot_border_right - res_x * 0.23 + heading_extents[3], plot_border_bottom - res_y * 0.02)
    show_text(cr, text_min)

    move_to(cr, plot_border_right - text_extents(cr, text_max)[3] - text_extents(cr, text_max)[1], plot_border_bottom - res_y * 0.02)
    show_text(cr, text_max)

    ## EXPORT SVG AND PNG
    #--------------------

    write_to_png(c, "wireplot.png")

    finish(c)

end

export wireplot

end # module
