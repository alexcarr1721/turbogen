

using HDF5
using DelimitedFiles
using Plots
gr()

# outputfile = ""
# outputname = "field"
# outputfile = string(outputname, string(1), ".h5")
outputfile = "turbo1.h5"
num_datasets = 500
x1 = h5read(outputfile, "x")
x2 = h5read(outputfile, "y")
x3 = h5read(outputfile, "z")
ux1_x1x2 = zeros(typeof(x1[1]), size(x1,1), size(x2,1))
ux1_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
ux1_x1x3 = zeros(typeof(x1[1]), size(x1,1), size(x3,1))
ux2_x1x2 = zeros(typeof(x1[1]), size(x1,1), size(x2,1))
ux2_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
ux2_x1x3 = zeros(typeof(x1[1]), size(x1,1), size(x3,1))
ux3_x1x2 = zeros(typeof(x1[1]), size(x1,1), size(x2,1))
ux3_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
ux3_x1x3 = zeros(typeof(x1[1]), size(x1,1), size(x3,1))
T_x1x2 = zeros(typeof(x1[1]), size(x1,1), size(x2,1))
T_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
T_x1x3 = zeros(typeof(x1[1]), size(x1,1), size(x3,1))

T_x1x2 = h5read(outputfile, "T", (255,:,:))
ux1_x1x2 = h5read(outputfile, "ux", (255,:,:))

T_plot = heatmap(x1, x2, T_x1x2)
display(T_plot)

ux_plot = heatmap(x1, x2, ux1_x1x2)
display(ux_plot)
