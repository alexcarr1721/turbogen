using HDF5
using LaTeXStrings
using DelimitedFiles
using Plots
gr()

# Load data ####################################################################
filename = "cbc_statistics.h5"
x1  = h5read(filename, "grid/x1")
x2  = h5read(filename, "grid/x2")
u1  = h5read(filename, "fields/u1")
u2  = h5read(filename, "fields/u2")
L_f = h5read(filename, "cbc/L_f")
cbc_std = h5read(filename, "cbc/std")
u0 = h5read(filename, "cbc/u0")
f_model = h5read(filename, "correlations/f_model")
g_model = h5read(filename, "correlations/g_model")
Longitudinal_correlation = h5read(filename, "correlations/Longitudinal_correlation")
Lateral_correlation = h5read(filename, "correlations/Lateral_correlation")
################################################################################

try
    mkdir("fig")
catch
    # Do nothing
end

# Grid length for cbc data
x_m = 0.0254 # meters

# read cbc data
cbc_longitudinal = readdlm("cbc_longitudinal.txt")
cbc_lateral = readdlm("cbc_lateral.txt")

# Contour plots ################################################################
u1_contour_fig = heatmap(x1./x_m, x2./x_m, u1'./u0,
    xlabel=latexstring("x_1 x_m^{-1}"), ylabel=latexstring("x_2 x_m^{-1}"),
    xlims=(-5,5), ylims=(-5,5),
    colorbar_title=latexstring("u_1 \\overline{u}^{-1}"), right_margin=4Plots.mm,
    xguidefontsize=14, yguidefontsize=14, colorbar_titlefontsize=14)
display(u1_contour_fig)
savefig(u1_contour_fig,"fig/u1_cbc.pdf")

u2_contour_fig = heatmap(x1./x_m, x2./x_m, u2'./u0,
    xlabel=latexstring("x_1 x_m^{-1}"), ylabel=latexstring("x_2 x_m^{-1}"),
    xlims=(-5,5), ylims=(-5,5),
    colorbar_title=latexstring("u_2 \\overline{u}^{-1}"), right_margin=4Plots.mm,
    xguidefontsize=14, yguidefontsize=14, colorbar_titlefontsize=14)
display(u2_contour_fig)
savefig(u2_contour_fig,"fig/u2_cbc.pdf")
################################################################################

# Grid length for cbc data
x_m = 0.0508 # meters

# Correlations #################################################################
fg_x1_plot = plot((x1.-x1[1])./x_m, Longitudinal_correlation,
    xlabel=latexstring("\\Delta x_1 x_m^{-1}"),
    ylabel=latexstring("f(\\Delta x_1 x_m^{-1}), \\;\\; g(\\Delta x_1 x_m^{-1})"),
    label=latexstring("f(\\Delta x_1 x_m^{-1})"), lw = 2, xlims=(0,6),
    ylims=(-0.25,1), linecolor=1,legendfontsize=12, xguidefontsize=14, yguidefontsize=14)
plot!(fg_x1_plot, (x1.-x1[1])./x_m, Lateral_correlation,
    label=latexstring("g(\\Delta x_1 x_m^{-1})"), lw = 2, linecolor=2)
plot!(fg_x1_plot, (x1.-x1[1])./x_m, f_model,
    label=latexstring("f(\\Delta x_1 x_m^{-1}) \\;\\; von \\; K\\'arm\\'an"),
    lw = 2, linestyle=:dash, linecolor=1)
plot!(fg_x1_plot, (x1.-x1[1])./x_m, g_model,
    label=latexstring("g(\\Delta x_1 x_m^{-1}) \\;\\; von \\; K\\'arm\\'an"),
    lw = 2, linestyle=:dash, linecolor=2)
scatter!(fg_x1_plot, cbc_longitudinal[2:end,1], cbc_longitudinal[2:end,2],
    label=latexstring("f(\\Delta x_1 x_m^{-1}) \\;\\; CBC \\; measured"),
    color=1)
scatter!(fg_x1_plot, cbc_lateral[2:end,1], cbc_lateral[2:end,2],
    label=latexstring("g(\\Delta x_1 x_m^{-1}) \\;\\; CBC \\; measured"),
    color=2)
display(fg_x1_plot)
savefig(fg_x1_plot, "fig/fg_cbc.pdf")
################################################################################
