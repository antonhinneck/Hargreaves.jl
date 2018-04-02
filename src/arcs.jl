
## DEFINING TYPE "ARCS"
#________________________________________

struct arcs

    arc_array::Array{Int64,2}
    arc_repo::Array{Int64,2}
    node_id_range::Tuple{Int,Int}
    arcs_df::DataFrame
    max_cons::Int64
    labels::Array{Symbol,1}
end

function arcs(arc_array::Array{Int64,2}, labels=[Symbol("")])

    node_id_range = extrema(arc_array)
    arc_repo = zeros(Int, node_id_range[2], node_id_range[2])
    arc_repo = store_arcs_to_repo(arc_array, arc_repo)

    arcs_df = from_repo_to_df(arc_repo)

    max_cons = get_max_cons(arcs_df)

    return arcs(arc_array, arc_repo, node_id_range, arcs_df, max_cons, labels)

end

function store_arcs_to_repo(arc_array::Array{Int64,2}, arc_repo)
    for i in 1:length(arc_array[:,1])
        arc_repo[arc_array[i,1], arc_array[i,2]] += 1
    end
    return arc_repo
end

get_max_cons(df::DataFrame) = maximum(df[:,3])
    
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
        for j in (i + 1):length(arc_repo[1,:])
            if arc_repo[i,j] > 0 || arc_repo[j,i] > 0
                push!(arcs_df, [i, j, arc_repo[i,j] + arc_repo[j,i], arc_repo[i,j], arc_repo[j,i]])
            end
        end
    end
    return arcs_df
end

