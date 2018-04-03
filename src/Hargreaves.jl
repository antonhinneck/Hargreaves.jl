__precompile__(true)
module Hargreaves

using Cairo, Colors, LightGraphs

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

function _get_color(value, max)
    (r_min, r_max) = (200, 112)
    (g_min, g_max) = (255, 30)
    (b_min, b_max) = (225, 112)
    r_out = (r_min - ((r_min - r_max) / max) * value) / 255
    g_out = (g_min - ((g_min - g_max) / max) * value) / 255
    b_out = (b_min - ((b_min - b_max) / max) * value) / 255
    return [
        [r_out,g_out,b_out], 
        [0.0,0.0,0.0], 
        [r_min / 255,g_min / 255,b_min / 255],
        [r_max / 255,g_max / 255,b_max / 255]
    ]
end


function wireplot(g::AbstractGraph, basefn = "wireplot"; 
    distmx=weights(g),
    res_x = 4096, res_y = 4096,
    plot_to_surface_ratio = 0.9,
    wireplot_bg_color=[1.0,1.0,1.0],
    wireplot_jitter=0,
    wireplot_node_diameter=6.0,
    wireplot_node_label_offset=16.0,
    wireplot_node_label_font_color=[0.0,0.0,0.0],
    wireplot_legend_font_color=[0.0,0.0,0.0],
    min_width=1,  # minimum normalized arc width
    max_width=100 # maximum normalized arc width
    )

    (actual_minwidth, actual_maxwidth) = extrema(distmx)
    actual_widthrange = actual_maxwidth - actual_minwidth
    @inline function _norm_width(w)
        if actual_widthrange == 0
            return 1
        else
            w_pct = (w - actual_minwidth) / actual_widthrange
            normrange = max_width - min_width
            if w_pct <= 0 || normrange <= 0
                return 1
            end
            return floor(Int, normrange * w_pct)
        end
    end
    
    max_cons = Δ(g)
    center = [res_x / 2,res_y / 2]

    c = CairoSVGSurface("$(basefn).svg", res_x, res_y)
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
    
    node_pos = _compute_node_positions(res_x, res_y, nv(g), plot_to_surface_ratio, 1)
    degs = degree(g)

    for (i, e) in enumerate(edges(g))
        (s, d) = (src(e), dst(e))
        # indeg_d = indegs[d]
        # deg_d = indeg_d + outdeg_d
        set_source_rgb(cr, _get_color(degs[d], max_cons)[1]...)
        set_line_width(cr, _norm_width(distmx[s, d]))
        curve_to(cr, node_pos[s, 1], node_pos[s, 2],
            center[1] + rand() * wireplot_jitter - rand() * wireplot_jitter,
            center[2] + rand() * wireplot_jitter - rand() * wireplot_jitter,
            node_pos[d, 1], node_pos[d, 2])
        stroke(cr)
    end
    set_line_width(cr, 1)

    ## PLOT NODES BASED ON COMPUTED POSITIONS
    #----------------------------------------

    for i in vertices(g)
        (pos_x, pos_y) = node_pos[i,:]
        
        set_source_rgb(cr, _get_color(degs[i], max_cons)[2]...)
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

    pattern_add_color_stop_rgb(gradient, 0, _get_color(0, max_cons)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, _get_color(max_cons, max_cons)[4]...)
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
    write_to_png(c, "$(basefn).png")
    finish(c)
        
end

export wireplot

end # module
