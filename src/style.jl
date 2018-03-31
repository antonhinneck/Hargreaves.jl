## VARIABLES
#-----------

wireplot_bg_color = [1.0,1.0,1.0]
wireplot_jitter = 0
wireplot_node_diameter = 6.0
wireplot_node_label_offset = 16.0
wireplot_node_label_font_color = [0.0,0.0,0.0]
wireplot_legend_font_color = [0.0,0.0,0.0]


## FUNCTIONS
#-----------

function get_color(value, max, template = Symbol("blossom"))

    if template == :blossom

       r_min = 200
       r_max = 112

       g_min = 255
       g_max = 30

       b_min = 225
       b_max = 112

       r_out = (r_min - ((r_min-r_max)/max)*value)/255
       g_out = (g_min - ((g_min-g_max)/max)*value)/255
       b_out = (b_min - ((b_min-b_max)/max)*value)/255

       return [[r_out,g_out,b_out], [0.0,0.0,0.0], [r_min/255,g_min/255,b_min/255], [r_max/255,g_max/255,b_max/255]]

    end
end
