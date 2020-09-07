# Post-processing code for input file

using HDF5
using Plots
using SpecialFunctions
using StatsBase
using Statistics
using Trapz
gr()

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

# Main Code

ux = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "ux", (40,:,:)) # Read slice
uy = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "uy", (40,:,:)) # Read slice
uz = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "uz", (40,:,:)) # Read slice
x = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "x") # Read slice
y = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "y") # Read slice
z = h5read("/home/ancarr/personal/tests/turb_generation/rgm_gaussian.h5",
    "z") # Read slice

uxu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "ux", (40,:,:)) # Read slice
uyu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "uy", (40,:,:)) # Read slice
uzu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "uz", (40,:,:)) # Read slice
xu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "x") # Read slice
yu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "y") # Read slice
zu = h5read("/home/ancarr/personal/tests/turb_generation/rgm.h5",
    "z") # Read slice

# Plots of turbulent field
uxplot = contour(x, y, ux, fill = true, linewidth = 0, xlabel = "x [cm]",
    label = "u_1")
yaxis!(uxplot, "y [cm]")
display(uxplot)
savefig(uxplot, "rpm.pdf")

uxplotu = contour(xu, yu, uxu, fill = true, linewidth = 0, xlabel = "x [cm]",
    label = "u_1")
yaxis!(uxplotu, "y [cm]")
display(uxplotu)
savefig(uxplotu, "rfm.pdf")


# Mean, Stddev
w = uweights(size(ux,1)*size(ux,2))
ux_mean = mean(ux,w)
ux_var = var(ux,w)

# Correlations
dx = x[2] - x[1]
M = 10/340
t_size = Int( floor(x[end]) ) + 1
ux_cor_stor = zeros(Float32, size(ux,1), t_size)
for i ∈ 1:size(ux,1)
    global τ
    ux_cor_stor[i,:], τ = autocorrelation(ux[i,:], x, x[end], ux_var)
end
ux_cor = dropmean(ux_cor_stor, dims=1)

# Mean, Stddev
wu = uweights(size(uxu,1)*size(uxu,2))
ux_meanu = mean(uxu,wu)
ux_varu = var(uxu,wu)

# Correlations
dxu = xu[2] - xu[1]
M = 10/340
t_sizeu = Int( floor(xu[end]) ) + 1
ux_cor_storu = zeros(Float32, size(uxu,1), t_sizeu)
for i ∈ 1:size(uxu,1)
    global τu
    ux_cor_storu[i,:], τu = autocorrelation(uxu[i,:], xu, xu[end], ux_varu)
end
ux_coru = dropmean(ux_cor_storu, dims=1)

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

τ[1] = 0.001
tm = τ./M
xgrid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
    100.0, 200.0, 300.0]
xshow = ["10⁻¹","","","","","","","","", "10⁰","","","","","","","","", "10¹",
    "","","","","","","","", "10²","",""]
# Plots
uxplot_cor = plot(tm,ux_cor,linewidth = 2, label = "Computation",
  linecolor = :steelblue, ylabel="R₁₁ (x₀M⁻¹=42, 0, 0; ΔxM⁻¹)", xaxis = :log10,
  xlabel=" ΔxM⁻¹ [cm] ",ylim = (0,1), xlim=(0.1,300), legend = :topright,
  xticks = (xgrid, xshow))
scatter!(uxplot_cor, xR11, R11, linewidth = 2, label = "CBC Measurement",
  color = :black)
display(uxplot_cor)
savefig(uxplot_cor,"Autocorrelationlog.pdf")

τ[1] = 0.0
uxplot_cor1 = plot(tm,ux_cor,linewidth = 2, label = "Computation",
  linecolor = :steelblue, ylabel="R₁₁ (x₀M⁻¹ = 42, 0, 0; ΔxM⁻¹)",
  xlabel=" ΔxM⁻¹ [cm] ",ylim = (0,1), xlim=(0,20), legend = :topright)
scatter!(uxplot_cor1, xR11, R11, linewidth = 2, label = "CBC Measurement",
  color = :black)
display(uxplot_cor1)
savefig(uxplot_cor1,"Autocorrelation.pdf")

# τu[1] = 0.001
# tmu = τu./M
# # Plots
# plot!(uxplot_cor, tm,ux_cor, linewidth = 2, label = "Uncorrelated phase",
#   color = :red)
# display(uxplot_cor)
# savefig(uxplot_cor,"Autocorrelationlog.pdf")
#
# τ[1] = 0.0
# plot!(uxplot_cor1, tm,ux_cor, linewidth = 2, label = "Uncorrelated phase",
#   color = :red)
# display(uxplot_cor1)
# savefig(uxplot_cor1,"Autocorrelation.pdf")
