

using HDF5
using DelimitedFiles
using Plots
gr()

function forward_diff1(
    f::AbstractArray{<:Number,N} where N,
    x::AbstractArray{<:Number,N} where N;
    bc="lastpoint"
    )

    df = zeros(typeof(f[1]), size(f,1))
    if ( bc == "lastpoint" )
        for i ∈ 1:size(x,1)-1
            df[i] = (f[i+1] - f[i])/(x[2] - x[1])
        end
        df[size(f,1)] = df[size(f,1)-1]
    elseif ( bc == "zero" )
        for i ∈ 1:size(x,1)-1
            df[i] = (f[i+1] - f[i])/(x[2] - x[1])
        end
        df[size(f,1)] = 0.0
    elseif ( bc == "periodic" )
        for i ∈ 1:size(x,1)-1
            df[i] = (f[i+1] - f[i])/(x[2] - x[1])
        end
        df[size(f,1)] = df[1]
    else
        for i ∈ 1:size(x,1)-1
            df[i] = (f[i+1] - f[i])/(x[2] - x[1])
        end
        df[size(f,1)] = df[size(f,1)-1]
    end

    return df
end

# outputfile = ""
# outputname = "field"
# outputfile = string(outputname, string(1), ".h5")
outputfile = "field1.h5"
num_datasets = 500
x1 = h5read(outputfile, "grid/x1")
x2 = h5read(outputfile, "grid/x2")
x3 = h5read(outputfile, "grid/x3")
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
rhoturb_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
rhoturbx_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
rhoturby_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))
rhoturbz_x2x3 = zeros(typeof(x1[1]), size(x2,1), size(x3,1))

T_x1x2 = h5read(outputfile, "turb/temperature", (127,:,:))
ux1_x1x2 = h5read(outputfile, "turb/ux1", (127,:,:))
rhoturb_x2x3 = h5read(outputfile, "turb/rho", (127,:,:))
rhoturby_x2x3 = h5read(outputfile, "turb/drhodx2", (127,:,:))

T_plot = heatmap(x1, x2, T_x1x2)
display(T_plot)

ux_plot = heatmap(x1, x2, ux1_x1x2)
display(ux_plot)

rhoturb_x2x3_plot = heatmap(x2, x3, rhoturb_x2x3)
display(rhoturb_x2x3)

rhoturby_x2x3_plot = heatmap(x2, x3, rhoturby_x2x3)
display(rhoturby_x2x3_plot)

# Compare derivatives?
drhodx2_forward = forward_diff1(rhoturb_x2x3[:,127],x2,bc="periodic")
drhodx2_fft     = rhoturby_x2x3[:,127]

L2norm = sqrt( sum((drhodx2_forward .- drhodx2_fft).^2) )

compare_plot = plot(x2, drhodx2_forward)
plot!(compare_plot, x2, drhodx2_fft)
display(compare_plot)
