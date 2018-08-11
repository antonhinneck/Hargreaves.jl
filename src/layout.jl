type layout

    #PLOT DIMENSIONS
    #The layout's components - being represented by rectangles - are defined by the left upper most point as well as their width and height respectively.

    canvas_dimensions::Array{Float64,1}
    plot_pane_def::Array{Float64, 1}
    legend_def::Array{Float64, 1}
    label_area_dimensions::Array{Float64, 1}

    #LAYOUT PROPERTIES
    show_legend::Bool
    show_labels::Bool

    function layout(canvas_dimensions::Array{Int64,1};
                    show_legend = true,
                    show_labels = true,
                    label_area_share = 0.05,
                    plot_pane_y_share = 0.7)

        #INTIT FIELDS
        plot_pane_def = Array{Float64, 1}(4)
        legend_def = Array{Float64, 1}(4)
        label_area_dimensions = Array{Float64, 1}(2)

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

        #Define components
        plot_pane_def = [canvas_dimensions[1] * label_area_share, canvas_dimensions[2] * label_area_share, canvas_dimensions[1] - canvas_dimensions[1] * 2 * label_area_share, canvas_dimensions[2] * plot_pane_y_share]
        legend_def = [canvas_dimensions[1] * label_area_share, plot_pane_def[4] + 2 *  label_area_share * canvas_dimensions[2], canvas_dimensions[1] - canvas_dimensions[1] * 2 * label_area_share, canvas_dimensions[2] - (plot_pane_def[4] + 2 * label_area_share * canvas_dimensions[2])]
        label_area_dimensions = [plot_pane_def[1], plot_pane_def[2]]
        new(canvas_dimensions, plot_pane_def, legend_def, label_area_dimensions, show_legend, show_labels)
    end
end
