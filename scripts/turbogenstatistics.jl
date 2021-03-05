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

ux1_x2x1_plot = heatmap(x₂, x₁, ux1_x2x1[:,:,1], xlabel="x₂ [m]",
    ylabel="x₁ [m]", colorbar_title="uₓ₁ [m/s]")
display(ux1_x2x1_plot)
savefig(ux1_x2x1_plot,"fig/ux1_x2x1_heatmap.pdf")

ux2_x2x1_plot = heatmap(x₂, x₁, ux2_x2x1[:,:,1], xlabel="x₂ [m]",
    ylabel="x₁ [m]", colorbar_title="uₓ₂ [m/s]")
display(ux2_x2x1_plot)
savefig(ux2_x2x1_plot,"fig/ux2_x2x1_heatmap.pdf")

ux3_x2x1_plot = heatmap(x₂, x₁, ux3_x2x1[:,:,1], xlabel="x₂ [m]",
    ylabel="x₁ [m]", colorbar_title="uₓ₃ [m/s]")
display(ux3_x2x1_plot)
savefig(ux3_x2x1_plot,"fig/ux3_x2x1_heatmap.pdf")

T_x2x1_plot = heatmap(x₂, x₁, T_x2x1[:,:,1], xlabel="x₂ [m]",
    ylabel="x₁ [m]", colorbar_title="T [K]")
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
################################################################################

# Sampling Frequency ###########################################################
fs = size(τ,1)/(τ[size(τ,1)] - τ[1])
################################################################################

# Correlations #################################################################

# Energy spectrum ##############################################################
