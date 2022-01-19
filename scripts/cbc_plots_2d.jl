using HDF5
using LaTeXStrings
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

# Contour plots ################################################################
u1_contour_fig = heatmap(x1./L_f, x2./L_f, u1'./u0,
    xlabel=latexstring("x_1 L_f^{-1}"), ylabel=latexstring("x_2 L_f^{-1}"),
    xlims=(-5,5), ylims=(-5,5), colorbar_title=latexstring("u_1 U_0^{-1}"))
display(u1_contour_fig)
savefig(u1_contour_fig,"fig/u1_cbc.pdf")

u2_contour_fig = heatmap(x1./L_f, x2./L_f, u2'./u0,
    xlabel=latexstring("x_1 L_f^{-1}"), ylabel=latexstring("x_2 L_f^{-1}"),
    xlims=(-5,5), ylims=(-5,5), colorbar_title=latexstring("u_2 U_0^{-1}"))
display(u2_contour_fig)
savefig(u2_contour_fig,"fig/u2_cbc.pdf")
################################################################################

# Correlations #################################################################
fg_x1_plot = plot(x1.-x1[1], Longitudinal_correlation, xlabel="r [m]",
    label="f(r)", lw = 2, xlims=(0,0.3), linecolor=1, ylabel="f(r), g(r)")
plot!(x1.-x1[1], Lateral_correlation, label="g(r)", lw = 2, linecolor=2)
plot!(x1.-x1[1], f_model, label="f(r) model", lw = 2, linestyle=:dash, linecolor=1)
plot!(x1.-x1[1], g_model, label="g(r) model", lw = 2, linestyle=:dash, linecolor=2)
display(fg_x1_plot)
savefig(fg_x1_plot, "fig/fg_x1.pdf")
################################################################################
