type layout

#LAYOUT PROPERTIES
show_legend
show_labels

#PLOT DIMENSIONS
canvas_dimensions
plot_pane_dimensions
label_area_dimensions
legend_dimensions

    function layout(; show_legend = false, show_labels = true)

        #INIT ARRAYS
        canvas_dimensions = Array{Int64, 1}(2)
        plot_pane_dimensions = Array{Int64, 1}(2)
        label_area_dimensions = Array{Int64, 1}(2)
        legend_dimensions = Array{Int64, 1}(2)

        for i in 1:2
            canvas_dimensions[i] = 0
            plot_pane_dimensions[i] = 0
            label_area_dimensions[i] = 0
            legend_dimensions[i] = 0
        end



    end

end
