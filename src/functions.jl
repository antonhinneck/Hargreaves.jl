function compute_node_positions(_res_x, _res_y, _node_count, ratio = 0.8, switch = 1)

    _center = [_res_x/2,_res_y/2]
    _degrees = 2*pi
    _degree_sections = _degrees/_node_count
    _main_radius = _res_x/2*ratio
    #println(string(_degree_sections," - ",_main_radius))

    _positions = Array{Float64}(_node_count,2)

    for i in 1:_node_count

        _positions[i,1] = cos(_degree_sections*(i-1))*_main_radius + _center[1]
        _positions[i,2] = sin(_degree_sections*(i-1))*_main_radius + _center[2]

    end
    return _positions
end
