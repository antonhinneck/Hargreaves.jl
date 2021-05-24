using LightGraphs
# using Hargreaves
include("C:/Users/Anton Hinneck/Desktop/Libya/Hargreaves.jl/src/Hargreaves.jl")
using .Hargreaves
# dir = @__DIR__
# cd(dir)

# Amount of vertices
# node_count = 31

# Upper triangular matrix with random values
# weights = Array{Int64, 2}(node_count,node_count)

include("C:/Users/Anton Hinneck/juliaPackages/GitHub/PowerGrids.jl/src/PowerGrids.jl")
using .PowerGrids
dir = @__DIR__
cd(dir)

set_csv_path("C:/Users/Anton Hinneck/Documents/Git/pglib2csv/pglib/2020-08-21.19-54-30-275/csv")
PowerGrids.csv_cases(verbose = true)
#PowerGrids.select_csv_case(48)
PowerGrids.select_csv_case(3)
case = PowerGrids.loadCase() # 118 Bus ieee
g = toGraph(case)

_weights = zeros(length(vertices(g)), length(vertices(g)))
for i in vertices(g)

    for j in vertices(g)

        _weights[j,i] = 0

        if j < i

            _weights[j,i] = 0
            _weights[j,i] = 1#round(Int, rand() * 20 + 1)

        else

            _weights[j,i] = 0

        end
    end
end

# Some random relatively high weight
# weights[2,15] = 30

# Define graph
#g = CompleteGraph(node_count)
cd(dir)
wireplot(g, "wireplot_1", distmx = _weights, wireplot_edge_color_scheme = :purple, static_widths = true, static_colors = true)
wireplot(g, "wireplot_2", distmx = _weights, wireplot_edge_color_scheme = :purple, color_scale = :quadratic, width_scale = :quadratic, max_width = 8)
