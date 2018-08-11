type layout

    #PLOT DIMENSIONS
    #The layout's components - being represented by rectangles - are defined by the left upper most point as well as their width and height respectively.

    canvas_dimensions::Array{Int64,1}
    plot_pane_def::Array{Int64, 1}
    legend_def::Array{Int64, 1}
    label_area_dimensions::Array{Int64, 1}

    #LAYOUT PROPERTIES
    show_legend::Bool
    show_labels::Bool

    function layout(canvas_dimensions::Array{Int64,1}; show_legend = true, show_labels = true)

        #INTIT FIELDS
        plot_pane_def = Array{Int64, 1}(4)
        legend_def = Array{Int64, 1}(4)
        label_area_dimensions = Array{Int64, 1}(2)

        for i in 1:4
            if i < 3
                label_area_dimensions[i] = 0
            end
            plot_pane_def[i] = 0
            legend_def[i] = 0
        end

        canvas_dimensions = canvas_dimensions
        show_legend = show_legend
        show_labels = show_labels

        #Defin

        new(canvas_dimensions, plot_pane_def, legend_def, label_area_dimensions, show_legend, show_labels)
    end
end

test = layout([200,200], show_labels = false, show_legend = true)
