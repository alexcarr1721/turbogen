# Statistics of turbulent field

push!(LOAD_PATH, pwd())

using DelimitedFiles
using FFTW
using Interpolations
using HDF5
using StatsBase
using Distributions
using Statistics
using KernelDensity
using LinearAlgebra
using Plots
gr()

# Functions
function autocorrelation(
    f::AbstractArray{<:Number,N} where N, # Signal
    x::AbstractArray{<:Number,N} where N  # Domain
    )

    # x must be evenly spaced
    dx = x[2] - x[1]

    # Compute mean and variance from Statistics
    f_mean = Statistics.mean(f)
    f_var  = Statistics.var(f)

    # Compute autocorrelation

    # zero pad by length of array
    n = size(f,1)
    fnew = zeros(typeof(f[1]), n*2)
    fnew[1:n] = f .- f_mean

    # Allocate autocorrelation array
    B = zero(f)

    # Multiply by delay of signal and average
    for i ∈ 1:n
        j = i-1
        B[i] = (1.0/(abs(x[n] - x[1])) )*sum(
            fnew[1:n].*fnew[1+j:n+j].*dx )/f_var
    end

    return B # return the normalized autocorrelation
end

function spectrum(
    f::AbstractArray{<:Number,N} where N, # Signal for many datasets
    x::AbstractArray{<:Number,N} where N  # Domain
    )

    # Sampling frequency
    fs = size(x,1)/(x[size(x,1)] - x[1])
    S  = zero(f[:,1])
    m  = size(f,2)
    n  = size(f,1)

    for i ∈ 1:n
        S[i] = (n/(fs*m))*sum((abs.(f[i,:])).^2)
    end

    return S # return the spectrum
end

function d2fdx2_6(
    f::AbstractArray{<:Number,N} where N, # Function
    x::AbstractArray{<:Number,N} where N  # Domain
    )

    # Automatic zero boundary condition
    dfdx = zero(f)
    for i ∈ 1:size(f,1)-6
        h = x[i+1] - x[i]
        dfdx[i] = (812*f[i]-3132*f[i+1]+5265*f[i+2]-5080*f[i+3]
            +2970*f[i+4]-972*f[i+5]+137*f[i+6])/(180*1.0*h^2)
    end

    return dfdx
end

function integrate_to_zero(
    f::AbstractArray{<:Number,N} where N, # Function
    x::AbstractArray{<:Number,N} where N  # Domain
    )

    index_zero = 0
    for i ∈ 1:size(f,1)-1
        if ( (f[i] > 0) && (f[i+1] <= 0) )
            index_zero = i + 1
            break
        end
    end

    integral = trapz(f[1:index_zero], x[1:index_zero])
    L = integral[index_zero]

    return L
end

function trapz(
    f::AbstractArray{<:Number,N} where N, # Function
    x::AbstractArray{<:Number,N} where N  # Domain
    )

    # Initialize Integral array
    Integral = zeros(Float64,size(f,1))

    # Compute cumulative integral along x
    for i = 1:size(f,1)
        if i == 1
            Integral[i] = 0.0
        else
            Integral[i] = (f[i] + f[i-1])*(x[i] - x[i-1])/2.0  + Integral[i-1]
        end
    end

    return Integral
end

function fourier_space(
    x::AbstractArray{<:Number,N} where N
    )

    f = zero(x)
    N = size(x,1)
    fs = N/( abs(x[N] - x[1]) )
    if mod(N,2) == 0
        f = (-N/2:N/2-1)*(fs/N)
    else
        f = (-(N-1)/2:(N-1)/2)*(fs/N)
    end

    f = fftshift(f)
    return f
end

# List of things to do
#   1) Plot Contour plots (check)
#   2) Compute correlations in each direction
#   3) Compute scalar statistics (check)
#   4) Compute energy spectrum
#   5) Compute sampling frequency (check)

try
    mkdir("fig")
catch
    println("Directory fig/ already exists.")
end

# Load data
filename = "TurbogenDataCase1.h5"
ux1_x2x3 = h5read(filename, "contours_x1_0/ux1")
ux2_x2x3 = h5read(filename, "contours_x1_0/ux2")
ux3_x2x3 = h5read(filename, "contours_x1_0/ux3")
T_x2x3 = h5read(filename, "contours_x1_0/T")
ux1_x2x1 = h5read(filename, "contours_x3_0/ux1")
ux2_x2x1 = h5read(filename, "contours_x3_0/ux2")
ux3_x2x1 = h5read(filename, "contours_x3_0/ux3")
T_x2x1 = h5read(filename, "contours_x3_0/T")
ux1_x1_centerline = h5read(filename, "centerline_x1/ux1")
ux2_x1_centerline = h5read(filename, "centerline_x1/ux2")
ux3_x1_centerline = h5read(filename, "centerline_x1/ux3")
T_x1_centerline = h5read(filename, "centerline_x1/T")
ux1_x2_centerline = h5read(filename, "centerline_x2/ux1")
ux2_x2_centerline = h5read(filename, "centerline_x2/ux2")
ux3_x2_centerline = h5read(filename, "centerline_x2/ux3")
T_x2_centerline = h5read(filename, "centerline_x2/T")
ux1_x3_centerline = h5read(filename, "centerline_x3/ux1")
ux2_x3_centerline = h5read(filename, "centerline_x3/ux2")
ux3_x3_centerline = h5read(filename, "centerline_x3/ux3")
T_x3_centerline = h5read(filename, "centerline_x3/T")
τ   = h5read(filename, "grid/tau")
x₁  = h5read(filename, "grid/x1")
x₂  = h5read(filename, "grid/x2")
x₃  = h5read(filename, "grid/x3")

# Contour plots ################################################################
ux1_x2x3_plot = heatmap(x₂, x₃, ux1_x2x3[:,:,1], xlabel="x₂ [m]",
    ylabel="x₃ [m]", colorbar_title="uₓ₁ [m/s]")
display(ux1_x2x3_plot)
savefig(ux1_x2x3_plot,"fig/ux1_x2x3_heatmap.pdf")

ux2_x2x3_plot = heatmap(x₂, x₃, ux2_x2x3[:,:,1], xlabel="x₂ [m]",
    ylabel="x₃ [m]", colorbar_title="uₓ₂ [m/s]")
display(ux2_x2x3_plot)
savefig(ux2_x2x3_plot,"fig/ux2_x2x3_heatmap.pdf")

ux3_x2x3_plot = heatmap(x₂, x₃, ux3_x2x3[:,:,1], xlabel="x₂ [m]",
    ylabel="x₃ [m]", colorbar_title="uₓ₃ [m/s]")
display(ux3_x2x3_plot)
savefig(ux3_x2x3_plot,"fig/ux3_x2x3_heatmap.pdf")

T_x2x3_plot = heatmap(x₂, x₃, T_x2x3[:,:,1], xlabel="x₂ [m]",
    ylabel="x₃ [m]", colorbar_title="T [K]")
display(T_x2x3_plot)
savefig(T_x2x3_plot,"fig/T_x2x3_heatmap.pdf")

ux1_x2x1_plot = heatmap(x₁[:], x₂, ux1_x2x1[:,:,1], xlabel="x₁ [m]",
    ylabel="x₂ [m]", colorbar_title="uₓ₁ [m/s]", aspect_ratio=:equal)
display(ux1_x2x1_plot)
savefig(ux1_x2x1_plot,"fig/ux1_x2x1_heatmap.pdf")

ux2_x2x1_plot = heatmap(x₁, x₂, ux2_x2x1[:,:,1], xlabel="x₁ [m]",
    ylabel="x₂ [m]", colorbar_title="uₓ₂ [m/s]", aspect_ratio=:equal)
display(ux2_x2x1_plot)
savefig(ux2_x2x1_plot,"fig/ux2_x2x1_heatmap.pdf")

ux3_x2x1_plot = heatmap(x₁, x₂, ux3_x2x1[:,:,1], xlabel="x₁ [m]",
    ylabel="x₂ [m]", colorbar_title="uₓ₃ [m/s]", aspect_ratio=:equal)
display(ux3_x2x1_plot)
savefig(ux3_x2x1_plot,"fig/ux3_x2x1_heatmap.pdf")

T_x2x1_plot = heatmap(x₁, x₂, T_x2x1[:,:,1], xlabel="x₁ [m]",
    ylabel="x₂ [m]", colorbar_title="T [K]", aspect_ratio=:equal)
display(T_x2x1_plot)
savefig(T_x2x1_plot,"fig/T_x2x1_heatmap.pdf")
################################################################################

# Scalar Statistics ############################################################
ux1_mean = mean(ux1_x2x3)
ux2_mean = mean(ux2_x2x3)
ux3_mean = mean(ux3_x2x3)
T_mean = mean(T_x2x3)
ux1_std = moment(ux1_x2x3, 2)
ux2_std = moment(ux2_x2x3, 2)
ux3_std = moment(ux3_x2x3, 2)
T_std = moment(T_x2x3, 2)
ux1_skw = moment(ux1_x2x3, 3)
ux2_skw = moment(ux2_x2x3, 3)
ux3_skw = moment(ux3_x2x3, 3)
T_skw = moment(T_x2x3, 3)
ux1_kur = moment(ux1_x2x3, 4)
ux2_kur = moment(ux2_x2x3, 4)
ux3_kur = moment(ux3_x2x3, 4)
T_kur = moment(T_x2x3, 4)

# Velocity trace 1/3 <uᵢuᵢ>
vtrace = (1/3)*(ux1_std^2 + ux2_std^2 + ux3_std^2)
################################################################################

# Sampling Frequency ###########################################################
fs = size(τ,1)/(τ[size(τ,1)] - τ[1])
################################################################################

# Correlations #################################################################
ux1_autocor_x1 = zero(ux1_x1_centerline)
ux2_autocor_x1 = zero(ux2_x1_centerline)
ux3_autocor_x1 = zero(ux3_x1_centerline)
T_autocor_x1 = zero(T_x1_centerline)
ux1_autocor_x2 = zero(ux1_x2_centerline)
ux2_autocor_x2 = zero(ux2_x2_centerline)
ux3_autocor_x2 = zero(ux3_x2_centerline)
T_autocor_x2 = zero(T_x2_centerline)
ux1_autocor_x3 = zero(ux1_x3_centerline)
ux2_autocor_x3 = zero(ux2_x3_centerline)
ux3_autocor_x3 = zero(ux3_x3_centerline)
T_autocor_x3 = zero(T_x3_centerline)

for i ∈ 1:size(ux1_x1_centerline,2)
    ux1_autocor_x1[:,i] = autocorrelation(ux1_x1_centerline[:,i], x₁)
    ux2_autocor_x1[:,i] = autocorrelation(ux2_x1_centerline[:,i], x₁)
    ux3_autocor_x1[:,i] = autocorrelation(ux3_x1_centerline[:,i], x₁)
    T_autocor_x1[:,i] = autocorrelation(T_x1_centerline[:,i], x₁)
end
for i ∈ 1:size(ux1_x2_centerline,2)
    ux1_autocor_x2[:,i] = autocorrelation(ux1_x2_centerline[:,i], x₂)
    ux2_autocor_x2[:,i] = autocorrelation(ux2_x2_centerline[:,i], x₂)
    ux3_autocor_x2[:,i] = autocorrelation(ux3_x2_centerline[:,i], x₂)
    T_autocor_x2[:,i] = autocorrelation(T_x2_centerline[:,i], x₂)
end
for i ∈ 1:size(ux1_x3_centerline,2)
    ux1_autocor_x3[:,i] = autocorrelation(ux1_x3_centerline[:,i], x₃)
    ux2_autocor_x3[:,i] = autocorrelation(ux2_x3_centerline[:,i], x₃)
    ux3_autocor_x3[:,i] = autocorrelation(ux3_x3_centerline[:,i], x₃)
    T_autocor_x3[:,i] = autocorrelation(T_x3_centerline[:,i], x₃)
end

ux1_autocor_avg_x1 = zero(ux1_autocor_x1[:,1])
ux2_autocor_avg_x1 = zero(ux2_autocor_x1[:,1])
ux3_autocor_avg_x1 = zero(ux3_autocor_x1[:,1])
T_autocor_avg_x1 = zero(T_autocor_x1[:,1])
ux1_autocor_avg_x2 = zero(ux1_autocor_x2[:,1])
ux2_autocor_avg_x2 = zero(ux2_autocor_x2[:,1])
ux3_autocor_avg_x2 = zero(ux3_autocor_x2[:,1])
T_autocor_avg_x2 = zero(T_autocor_x2[:,1])
ux1_autocor_avg_x3 = zero(ux1_autocor_x3[:,1])
ux2_autocor_avg_x3 = zero(ux2_autocor_x3[:,1])
ux3_autocor_avg_x3 = zero(ux3_autocor_x3[:,1])
T_autocor_avg_x3 = zero(T_autocor_x3[:,1])

for i ∈ 1:size(ux1_autocor_x1,1)
    ux1_autocor_avg_x1[i] = sum(ux1_autocor_x1[i,:])/size(ux1_x1_centerline,2)
    ux2_autocor_avg_x1[i] = sum(ux2_autocor_x1[i,:])/size(ux2_x1_centerline,2)
    ux3_autocor_avg_x1[i] = sum(ux3_autocor_x1[i,:])/size(ux3_x1_centerline,2)
    T_autocor_avg_x1[i] = sum(T_autocor_x1[i,:])/size(T_x1_centerline,2)
end
for i ∈ 1:size(ux1_autocor_x2,1)
    ux1_autocor_avg_x2[i] = sum(ux1_autocor_x2[i,:])/size(ux1_x2_centerline,2)
    ux2_autocor_avg_x2[i] = sum(ux2_autocor_x2[i,:])/size(ux2_x2_centerline,2)
    ux3_autocor_avg_x2[i] = sum(ux3_autocor_x2[i,:])/size(ux3_x2_centerline,2)
    T_autocor_avg_x2[i] = sum(T_autocor_x2[i,:])/size(T_x2_centerline,2)
end
for i ∈ 1:size(ux1_autocor_x3,1)
    ux1_autocor_avg_x3[i] = sum(ux1_autocor_x3[i,:])/size(ux1_x3_centerline,2)
    ux2_autocor_avg_x3[i] = sum(ux2_autocor_x3[i,:])/size(ux2_x3_centerline,2)
    ux3_autocor_avg_x3[i] = sum(ux3_autocor_x3[i,:])/size(ux3_x3_centerline,2)
    T_autocor_avg_x3[i] = sum(T_autocor_x3[i,:])/size(T_x3_centerline,2)
end

# ux1_autocor_x1_plot = plot(x₁, ux1_autocor_avg_x1, xlabel="r₁ [m]",
#     ylabel="ρ₁₁(r₁)", label=:none)
# display(ux1_autocor_x1_plot)
# savefig(ux1_autocor_x1_plot, "fig/ux1_autocor_x1.pdf")
# ux2_autocor_x1_plot = plot(x₁, ux2_autocor_avg_x1, xlabel="r₁ [m]",
#     ylabel="ρ₁₁(r₁)", label=:none)
# display(ux2_autocor_x1_plot)
# savefig(ux2_autocor_x1_plot, "fig/ux2_autocor_x1.pdf")
# ux3_autocor_x1_plot = plot(x₁, ux3_autocor_avg_x1, xlabel="r₁ [m]",
#     ylabel="ρ₁₁(r₁)", label=:none)
# display(ux3_autocor_x1_plot)
# savefig(ux3_autocor_x1_plot, "fig/ux3_autocor_x1.pdf")
# T_autocor_x1_plot = plot(x₁, T_autocor_avg_x1, xlabel="r₁ [m]",
#     ylabel="ρ₁₁(r₁)", label=:none)
# display(T_autocor_x1_plot)
# savefig(T_autocor_x1_plot, "fig/T_autocor_x1.pdf")

fg_x1_plot = plot(x₁, ux1_autocor_avg_x1, xlabel="r [m]",
    ylabel="ρ₁₁(r₁)", label="f(r)")
plot!(x₁, ux2_autocor_avg_x1, label="g(r)")
display(fg_x1_plot)
savefig(fg_x1_plot, "fig/fg_x1.pdf")

# Taylor microscale ############################################################
fpp = d2fdx2_6(ux1_autocor_avg_x1, x₁)
gpp = d2fdx2_6(ux2_autocor_avg_x1, x₁)

λ = sqrt(-1.0/fpp[1])
mean_square_vorticity = 15*vtrace/(λ^2)

Lp = integrate_to_zero(ux1_autocor_avg_x1, x₁)
Ln = integrate_to_zero(ux2_autocor_avg_x1, x₁)

# Write to file
io = open("turbulent_field.txt", "w")
write(io, "up^2     un^2    Lp    Ln    Taylor Microscale   Mean-Square Vorticity")
write(io,"\n")
write(io, string(ux1_std^2, "   ",ux2_std^2,"   ",Lp,"   ",
    Ln,"   ",λ,"   ",mean_square_vorticity))
write(io,"\n")
write(io,"\n")
write(io,"\n")
write(io, "Mean ux1: "); write(io, ux1_mean)
write(io,"\n")
write(io, "Mean ux2: "); write(io, ux2_mean)
write(io,"\n")
write(io, "Mean ux3: "); write(io, ux3_mean)
write(io,"\n")
write(io, "Std dev ux1: "); write(io, ux1_std)
write(io,"\n")
write(io, "Std dev ux2: "); write(io, ux2_std)
write(io,"\n")
write(io, "Std dev ux3: "); write(io, ux3_std)
write(io,"\n")
write(io, "1/3(ux1_var + ux2_var + ux3_var): "); write(io, vtrace)
close(io)
# Energy spectrum ##############################################################

S₁₁ = spectrum(ux1_x1_centerline,x₁)
S₂₂ = spectrum(ux2_x1_centerline,x₁)
S₃₃ = spectrum(ux3_x1_centerline,x₁)
f₁  = fourier_space(x₁)
f₂  = fourier_space(x₂)
f₃  = fourier_space(x₃)

S11_plot = plot(f₁[2:Int(size(f₁,1)/2)],S₁₁[2:Int(size(f₁,1)/2)],xscale=:log10,yscale=:log10,xlabel="f [Hz]",ylabel="S₁₁",
    label=:none)
display(S11_plot)
savefig(S11_plot,"fig/S11_plot.pdf")
S22_plot = plot(f₁[2:Int(size(f₁,1)/2)],S₂₂[2:Int(size(f₁,1)/2)],xscale=:log10,yscale=:log10,xlabel="f [Hz]",ylabel="S₂₂",
    label=:none)
display(S22_plot)
savefig(S22_plot,"fig/S22_plot.pdf")
