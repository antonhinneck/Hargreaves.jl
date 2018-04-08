using LightGraphs
using Hargreaves

cd("C:/Users/anton/.julia/v0.6/Hargreaves/test")

# Amount of vertices
node_count = 31

# Upper triangular matrix with random values
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

# Some random relatively high weight
weights[2,15] = 30

# Define graph
g = CompleteGraph(node_count)

wireplot(g, "wireplot_1", distmx = weights, wireplot_edge_color_scheme = :purple, static_widths = true, static_colors = true)
wireplot(g, "wireplot_2", distmx = weights, wireplot_edge_color_scheme = :purple, color_scale = :quadratic, width_scale = :quadratic, max_width = 8)
