using Hargreaves
using Base.Test

cd("C:/Users/anton/.julia/v0.6/Hargreaves/test")

node1 = trunc(Int64, rand(1200)*32+1)
node2 = trunc(Int64, rand(1200)*32+1)

connections = hcat(node1,node2)

wireplot(connections)
