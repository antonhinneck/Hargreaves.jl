module Hargreaves

using Cairo, Colors, DataFrames

include("style.jl")
include("functions.jl")
include("arcs.jl")

function wireplot(connections::Array{Int64,2}, res_x = 4096, res_y = 4096, plot_to_surface_ratio = 0.9, plot_name = "wireplot.svg")

    data = arcs(connections)

    center = [res_x/2,res_y/2]
    border_left = center[1]
    border_right = 0.0
    border_top = 0.0
    border_bottom = 0.0

    c = CairoSVGSurface("wireplot.svg", res_x, res_y)
    cr = CairoContext(c)

    select_font_face(cr, "Times", 1,1)
    set_font_size(cr, 56.0)

    ## DRAW BACKGROUND
    #-----------------

    set_source_rgb(cr, wireplot_bg_color...)
    rectangle(cr, 0.0,0.0,res_x,res_y)
    fill(cr)

    ## COMPUTE NODE POSITIONS
    #------------------------

    node_positions = compute_node_positions(res_x,res_y,length(data.arc_repo[:,1]), plot_to_surface_ratio,1)

    ## PLOT NODE CONNECTIONS
    #-----------------------

    for i in 1:length(data.arcs_df[:,1])

        if data.arcs_df[i,1] > 0 || data.arcs_df[i,2] > 0

            set_source_rgb(cr, get_color(data.arcs_df[i,3], data.max_cons)[1]...)

            curve_to(cr, node_positions[data.arcs_df[i,1], 1],
                         node_positions[data.arcs_df[i,1], 2],
                         center[1] + rand()*wireplot_jitter - rand()*wireplot_jitter,
                         center[2] + rand()*wireplot_jitter - rand()*wireplot_jitter,
                         node_positions[data.arcs_df[i,2], 1],
                         node_positions[data.arcs_df[i,2], 2])
            stroke(cr)

        end

    end

    ## PLOT NODES BASED ON COMPUTED POSITIONS
    #----------------------------------------

    for i in 1:length(data.arc_repo[:,1])

        set_source_rgb(cr, get_color(data.arcs_df[1,3], data.max_cons)[2]...)
        arc(cr,node_positions[i,1], node_positions[i,2], wireplot_node_diameter, wireplot_node_diameter, 2*wireplot_node_diameter*pi)
        fill(cr)

    end

    ## PLOT NODES' LABELS
    #--------------------

    for i in 1:length(data.arc_repo[:,1])

        text = string(i)
        label_extents = text_extents(cr, text)

        set_source_rgb(cr, wireplot_node_label_font_color...)

        if round(Int,node_positions[i,1]) == center[1] && node_positions[i,2] > center[2]

            move_to(cr,node_positions[i,1] - label_extents[3]/2 - 1, node_positions[i,2] + wireplot_node_diameter/2 + label_extents[4] + wireplot_node_label_offset)
            show_text(cr, string(i))
            border_bottom = node_positions[i,2] + wireplot_node_diameter/2 + label_extents[4] + wireplot_node_label_offset

        elseif round(Int,node_positions[i,1]) == center[1] && node_positions[i,2] < center[2]

            move_to(cr,node_positions[i,1] - label_extents[3]/2 - 1, node_positions[i,2] - wireplot_node_diameter/2 - wireplot_node_label_offset)
            show_text(cr, string(i))
            border_top = node_positions[i,2] - wireplot_node_diameter/2 - wireplot_node_label_offset - text_extents(cr, string(i))[4]

        elseif node_positions[i,1] < center[1]

            move_to(cr,node_positions[i,1] - wireplot_node_diameter/2 - label_extents[3] - wireplot_node_label_offset, node_positions[i,2] + label_extents[4]/2 + sin((i*2*pi)/length(node_positions[:,1]))*(wireplot_node_label_offset + wireplot_node_diameter/2 + label_extents[4]/2))
            show_text(cr, string(i))

            if (node_positions[i,1] - wireplot_node_diameter/2 - label_extents[3] - wireplot_node_label_offset) < border_left

                border_left = node_positions[i,1] - wireplot_node_diameter/2 - label_extents[3] - wireplot_node_label_offset

            end

        elseif node_positions[i,1] > center[1]

            move_to(cr,node_positions[i,1] + wireplot_node_diameter/2 + wireplot_node_diameter/2 + wireplot_node_label_offset, node_positions[i,2] + label_extents[4]/2 + sin((i*2*pi)/length(node_positions[:,1]))*(wireplot_node_label_offset + wireplot_node_diameter/2 + label_extents[4]/2))
            show_text(cr, string(i))

            if (node_positions[i,1] + wireplot_node_diameter/2 + wireplot_node_diameter/2 + wireplot_node_label_offset + label_extents[3]) > border_right

                border_right = node_positions[i,1] + wireplot_node_diameter/2 + wireplot_node_diameter/2 + wireplot_node_label_offset + label_extents[3]

            end

        end
    end

    fill(cr)

    ## PLOT LEGEND
    #-------------

        ## HEADING

    println(string(border_top))
    println(string(border_left))
    println(string(border_right))
    println(string(border_bottom))

    heading = "Transported pieces"
    heading_extents = text_extents(cr, heading)

    set_source_rgb(cr, wireplot_legend_font_color...)
    move_to(cr, border_right - res_x*0.25, border_bottom)
    show_text(cr, heading)
    line_to(cr, border_right, border_bottom)

    stroke(cr)

        ## COLOR DESCRIPTOR

    rectangle(cr, border_right - res_x*0.23 + heading_extents[3],  border_bottom - res_y*0.015,  res_x*0.23 - heading_extents[3], res_y*0.01)

    gradient = pattern_create_linear(border_right - res_x*0.23 + heading_extents[3], border_bottom,  border_right,  border_bottom)

    pattern_add_color_stop_rgb(gradient, 0, get_color(0, data.max_cons)[3]...)
    pattern_add_color_stop_rgb(gradient, 1, get_color(data.max_cons, data.max_cons)[4]...)

    set_source(cr, gradient)

    fill_preserve(cr)

    destroy(gradient)

        ## COLOR DESCRIPTOR FONT

    set_source_rgb(cr, wireplot_legend_font_color...)

    text_min = string("[",1,", ...")
    text_max = string(data.max_cons,"]")

    move_to(cr, border_right - res_x*0.23 + heading_extents[3], border_bottom - res_y*0.02)
    show_text(cr, text_min)

    move_to(cr, border_right - text_extents(cr, text_max)[3] - text_extents(cr, text_max)[1], border_bottom - res_y*0.02)
    show_text(cr, text_max)

    ## EXPORT SVG AND PNG
    #--------------------

    write_to_png(c,"wireplot.png")

    finish(c)

end

export wireplot

end # module
