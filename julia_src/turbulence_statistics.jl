

using HDF5
using DelimitedFiles
using Statistics
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

# Functions
dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)

function autocorrelation(f::AbstractArray{<:Number,N} where N,
    t::AbstractArray{<:Number,N} where N,
    extension::Number,
    variance::Number)

    if (size(t,1) == 1)
        τ = 0.0
    else
        dt = t[2] - t[1]
        τ  = 0.0:dt:extension*dt
    end
    τ = convert(Array{Float64}, τ)

    f_new = zeros(Float64, size(f,1) + size(τ,1) )
    f_new[1:size(f,1)] = f[1:size(f,1)]

    # Construct integrand
    integrand = zero(f)
    integral  = zeros(Float64, size(τ,1))
    for i ∈ 1:length(τ)
        j = i - 1
        for k ∈ 1:size(f,1)
            integrand[k] = (1.0/t[end])*(f_new[k])*(f_new[k+j])/variance
        end
        integral[i] = trapz(t, integrand)
    end

    return integral, τ
end

# outputfile = ""
# outputname = "field"
# outputfile = string(outputname, string(1), ".h5")
outputfile = "field3.h5"
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

# T_x1x2 = h5read(outputfile, "turb/temperature", (127,:,:))
# ux1_x1x2 = h5read(outputfile, "turb/ux1", (127,:,:))
# rhoturb_x2x3 = h5read(outputfile, "turb/rho", (127,:,:))
# rhoturby_x2x3 = h5read(outputfile, "turb/drhodx2", (127,:,:))

# T_plot = heatmap(x1, x2, T_x1x2)
# display(T_plot)
#
# ux_plot = heatmap(x1, x2, ux1_x1x2)
# display(ux_plot)
#
# rhoturb_x2x3_plot = heatmap(x2, x3, rhoturb_x2x3)
# display(rhoturb_x2x3)
#
# rhoturby_x2x3_plot = heatmap(x2, x3, rhoturby_x2x3)
# display(rhoturby_x2x3_plot)

# Compare derivatives?
# drhodx2_forward = forward_diff1(rhoturb_x2x3[:,127],x2,bc="periodic")
# drhodx2_fft     = rhoturby_x2x3[:,127]
#
# L2norm = sqrt( sum((drhodx2_forward .- drhodx2_fft).^2) )
#
# compare_plot = plot(x2, drhodx2_forward)
# plot!(compare_plot, x2, drhodx2_fft)
# display(compare_plot)


# Statistics

ux1 = h5read(outputfile, "turb/ux1")
ux2 = h5read(outputfile, "turb/ux2")
ux3 = h5read(outputfile, "turb/ux3")
T   = h5read(outputfile, "turb/temperature")

w = uweights(size(ux1,1)*size(ux1,2)*size(ux1,3))
ux1_mean = mean(ux1,w)
ux1_var = var(ux1,w)

t_size = Int( floor(x1[end]) ) + 1
ux_cor_stor = zeros(Float32, size(ux1,2)*size(ux1,3), t_size)
for k ∈ 1:size(ux1,3)
    for j ∈ 1:size(ux1,2)
        global τ
        i = j + (k-1)*size(ux1,2)
        ux_cor_stor[i,:], τ = autocorrelation(ux1[:,j,k], x1, x1[end], ux1_var)
    end
end
ux_cor = dropmean(ux_cor_stor, dims=1)

# CBC data
R11 = [0.93786,
0.89069,
0.83088,
0.76431,
0.72093,
0.58112,
0.53603,
0.40445,
0.09982,
0.05636
]
xR11 = [0.38281,
0.77099,
1.31153,
2.5523,
4.0725,
6.0723,
8.1184,
17.987,
123.482,
223.56
]

# τ[1] = 0.001
# xm = 5.08
# tm = τ./xm
# xgrid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0,
#     6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
#     100.0, 200.0, 300.0]
# xshow = ["10⁻¹","","","","","","","","", "10⁰","","","","","","","","", "10¹",
#     "","","","","","","","", "10²","",""]
# # Plots
# uxplot_cor = plot(tm,ux_cor./ux_cor[1],linewidth = 2, label = "Computation",
#   linecolor = :steelblue, ylabel="R₁₁ (x₀/xₘ=42, 0, 0; Δx/xₘ)", xaxis = :log10,
#   xlabel=" Δx/xₘ",ylim = (0,1), xlim=(0.1,300), legend = :topright,
#   xticks = (xgrid, xshow))
# scatter!(uxplot_cor, xR11, R11, linewidth = 2, label = "CBC Measurement",
#   color = :black)
# display(uxplot_cor)
# savefig(uxplot_cor,"Autocorrelationlog.pdf")
#
# τ[1] = 0.0
# uxplot_cor1 = plot(tm,ux_cor./ux_cor[1],linewidth = 2, label = "Computation",
#   linecolor = :steelblue, ylabel="R₁₁ (x₀/xₘ = 42, 0, 0; Δx/xₘ)",
#   xlabel=" Δx/xₘ",ylim = (0,1), xlim=(0,20), legend = :topright)
# scatter!(uxplot_cor1, xR11, R11, linewidth = 2, label = "CBC Measurement",
#   color = :black)
# display(uxplot_cor1)
# savefig(uxplot_cor1,"Autocorrelation.pdf")
#
# ux_new_plot = plot(τ,ux_cor./ux_cor[1],linewidth = 2, label = "Computation",
#   linecolor = :steelblue, ylabel="R₁₁ (x₀/xₘ = 42, 0, 0; Δx)",
#   xlabel=" Δx ",ylim = (0,1), xlim=(0,20), legend = :topright)
# display(ux_new_plot)
# savefig(ux_new_plot,"Autocor_x.pdf")

Int_length = trapz(τ, ux_cor)
