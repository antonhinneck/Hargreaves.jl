using LightGraphs
using Hargreaves

cd("C:/Users/USERNAME/.julia/v0.6/Hargreaves/test")

node_count = 31

weights = Array{Int64, 2}(node_count,node_count)

for i in 1:node_count

    for j in 1:node_count

        weights[j,i] = 0

        if j < i

            weights[j,i] = 0
            weights[j,i] = round(Int, rand() * 20 + 1)

        else

            weights[j,i] = 0

        end
    end
end

weights[2,15] = 94

g = CompleteGraph(node_count)

wireplot(g, "wireplot_1", distmx = weights, wireplot_edge_color_scheme = :purple, max_width = 6)
wireplot(g, "wireplot_2", distmx = weights, wireplot_edge_color_scheme = :purple, static_colors = true, static_widths = true)
