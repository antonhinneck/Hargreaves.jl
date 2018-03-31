
## DEFINING TYPE "ARCS"
#________________________________________

type arcs

    arc_array::Array{Int64,2}
    arc_repo::Array{Int64,2}
    node_id_range::Array{Int64,1}
    arcs_df::DataFrame
    max_cons::Int64
    labels::Array{Symbol,1}

    function arcs(arc_array::Array{Int64,2}, labels = [Symbol("")])

        node_id_range = get_node_id_range(arc_array)
        arc_repo = init_repo(node_id_range)
        arc_repo = store_arcs_to_repo(arc_array, arc_repo)

        arcs_df = from_repo_to_df(arc_repo)

        max_cons = get_max_cons(arcs_df)

        new(arc_array, arc_repo, node_id_range, arcs_df, max_cons, labels)

    end

    function store_arcs_to_repo(arc_array::Array{Int64,2}, arc_repo)

        for i in 1:length(arc_array[:,1])

            arc_repo[arc_array[i,1], arc_array[i,2]] += 1

        end

        return arc_repo
    end

    function get_max_cons(df::DataFrame)

        max = 0

        for(i, value) in enumerate(df[:,3])

            if value > max

                max = value

            end
        end

        return max
    end

    function init_repo(node_id_range::Array{Int64,1})

        arc_repo = Array{Int64,2}(node_id_range[2], node_id_range[2])

        for i in 1:node_id_range[2]

            for j in 1:node_id_range[2]

                arc_repo[i,j] = 0

            end
        end

        return arc_repo

    end

    function get_node_id_range(arc_array::Array{Int64,2})

        min = arc_array[1,1]
        max = arc_array[1,1]

        for i in 1:length(arc_array[:,1])

            if arc_array[i,1] > max

                max = arc_array[i,1]

            end
            if arc_array[i,1] < min

                min = arc_array[i,1]

            end
            if arc_array[i,2] > max

                max = arc_array[i,2]

            end
            if arc_array[i,2] < min

                min = arc_array[i,2]

            end
        end

        return [min,max]
    end

    function from_repo_to_df(arc_repo)

        column_eltypes = Vector{DataType}()
        column_eltypes = [Int64, Int64, Int64, Int64, Int64]

        cnames = Vector{Symbol}()
        cnames = [Symbol("node1"), Symbol("node2"), Symbol("cons"),Symbol("to_node2"), Symbol("from:_node2")]

        ispda = Vector{Bool}()
        ispda = [false, false,false,false,false]
        ncols = 5

        arcs_df = DataFrame(column_eltypes, cnames, ispda, ncols)
        deleterows!(arcs_df, 1:5)

        for i in 1:length(arc_repo[:,1])

            for j in (i+1):length(arc_repo[1,:])

                if arc_repo[i,j] > 0 || arc_repo[j,i] > 0

                    push!(arcs_df, [i, j, arc_repo[i,j] + arc_repo[j,i], arc_repo[i,j], arc_repo[j,i]])

                end
            end
        end

        return arcs_df
    end
end
