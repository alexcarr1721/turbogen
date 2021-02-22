

using HDF5
using DelimitedFiles
using FFTW
using Statistics
using Trapz
using Plots
gr()

outputfile      = ""
outputname      = "field"
num_datasets    = 500
outputfile      = string(outputname, string(1), ".h5")
x               = h5read(outputfile, "x")
y               = h5read(outputfile, "y")
z               = h5read(outputfile, "z")
ux_xy           = zeros(typeof(x[1]), size(x,1), size(y,1), num_datasets)
ux_yz           = zeros(typeof(x[1]), size(y,1), size(z,1), num_datasets)
ux_xz           = zeros(typeof(x[1]), size(x,1), size(z,1), num_datasets)
uy_xy           = zeros(typeof(x[1]), size(x,1), size(y,1), num_datasets)
uy_yz           = zeros(typeof(x[1]), size(y,1), size(z,1), num_datasets)
uy_xz           = zeros(typeof(x[1]), size(x,1), size(z,1), num_datasets)
uz_xy           = zeros(typeof(x[1]), size(x,1), size(y,1), num_datasets)
uz_yz           = zeros(typeof(x[1]), size(y,1), size(z,1), num_datasets)
uz_xz           = zeros(typeof(x[1]), size(x,1), size(z,1), num_datasets)
T_xy            = zeros(typeof(x[1]), size(x,1), size(y,1), num_datasets)
T_yz            = zeros(typeof(x[1]), size(y,1), size(z,1), num_datasets)
T_xz            = zeros(typeof(x[1]), size(x,1), size(z,1), num_datasets)


count = 0
error = 0
for i ∈ 1:num_datasets
    global outputfile, x, y, z, ux_xy, uy_xy, uz_xy, T_xy

    outputfile = string(outputname, string(i), ".h5")
    T_xy[:,:,i] = h5read(outputfile, "T", (Int(size(z,1)/2 - 1),:,:))
    ux_xy[:,:,i] = h5read(outputfile, "ux", (Int(size(z,1)/2 - 1),:,:))
    uy_xy[:,:,i] = h5read(outputfile, "uy", (Int(size(z,1)/2 - 1),:,:))
    uz_xy[:,:,i] = h5read(outputfile, "uz", (Int(size(z,1)/2 - 1),:,:))

end


# T_x1x2 = h5read(outputfile, "T", (127,:,:))
# ux1_x1x2 = h5read(outputfile, "ux", (127,:,:))

# Calculate statistics
# Process the data to compute spectra and correlations
M = size(ux_xy,3) # Number of datasets to average over is equal to the size of dim3
N = size(ux_xy,1) # Number of points per dataset
transverse = size(ux_xy,2)
U = 10    # [m/s] freestream velocity to perform Taylor hypothesis
t = x./U    # Time series based on Taylor's Hypothesis
dt = t[2] - t[1]
Sᵤᵤ = zeros(typeof(ux_xy[1,1,1]), N, transverse)
u1  = zeros(ComplexF32, N )
Rᵤᵤ = zeros(typeof(ux_xy[1,1,1]), N, transverse)
fₛ = N/t[end]
f = fₛ*(0:(N/2))/N
c₀ = 343 # m/s
k = 2*π*f/c₀

for m ∈ 1:M
    for j ∈ 1:transverse
        u1[:]   =   fft(ux_xy[:,j,m]./100) # Convert to m/s
        for n ∈ 1:N
            Sᵤᵤ[n,j]   = Sᵤᵤ[n,j] + ( 1/(fₛ*M*N) )*abs(u1[n]^2)
        end
        Rᵤᵤ[:,j] = real.(ifft(Sᵤᵤ[:,j])).*N
    end
end

# Parameters for plotting
σᵤ = 0.222 # [m/s] from CBC
Lᵤ = 0.024 # [m] from CBC

# Comte-Bellot Corrsin data
E_2in = readdlm("cbc_e_2.txt", ',', '\n',skipstart=2)
replace!(E_2in, 0.0=>NaN)

σ_u_2in = [22.2, 12.8, 8.95]
ϵ_2in = [4740, 633, 174]
η_2in = [0.029, 0.048, 0.066]
λ_2in = [0.484, 0.764, 1.02]
L_2in = [1.27, 1.88, 2.28]
L_f_2in = [2.40, 3.45, 4.90]
R_λ_2in = [71.6, 65.3, 60.7]

R11_508 = readdlm("cbc_r11_508.txt", ',', '\n',skipstart=2)
R11_254 = readdlm("cbc_r11_254.txt", ',', '\n',skipstart=2)

Suu_plot = plot(f[2:end].*Lᵤ./U,Sᵤᵤ[2:Int((N/2 + 1)),100].*f[2:end]./(σᵤ^2),
    label=:none,xscale=:log10,yscale=:log10,ylabel="fS₁₁/σᵤ²",
    xlabel="fLᵤ/U")#,xticks = [0.1, 1, 10], xlim=(0.1,20), ylim=(0.1,20))
display(Suu_plot)
savefig(Suu_plot, string("Suu_plot.pdf"))

Suu_plotd = plot((2*π.*f[2:end]./c₀).*(1/100),Sᵤᵤ[2:Int((N/2 + 1)),100].*(100^3),
    label=:none,xscale=:log10,yscale=:log10,ylabel="fS₁₁/σᵤ²",
    xlabel="fLᵤ/U")
display(Suu_plotd)
savefig(Suu_plotd, string("Suu_plotd.pdf"))

E11_plot = plot((2*π.*f[2:end]./c₀).*100, Sᵤᵤ[2:Int((N/2 + 1)),100].*(100^3),
    label=:none,xscale=:log10,yscale=:log10,ylabel="E₁₁ cm³/s⁻²",
    xlabel="k₁ cm⁻¹", lw=2)
plot!(E11_plot, E_2in[:,1], E_2in[:,2], label=:none, linecolor=:black)
display(E11_plot)
savefig(E11_plot, string("E11_plot.pdf"))

Mach = 5.08/100

t_disp = (x .- x[1])./Mach#./U
t_data = R11_508[:,1]#.*(M/U)

Ruu_plot = plot(t_disp[1:75], Rᵤᵤ[1:75,100]./Rᵤᵤ[1,100],
    label="Generated Field", xlabel="UΔtM⁻¹", xlim=(0,200), ylim=(0,1),
    ylabel="R₁₁ (x₀M⁻¹ = 42, 0, 0; ΔxM⁻¹)",legendfontsize=8, lw=2)
scatter!(Ruu_plot, t_data, R11_508[:,2],label="Comte-Bellot and Corrsin (1971)")
display(Ruu_plot)
savefig(Ruu_plot, string("Ruu_plot.pdf"))

Ruu_plot1 = plot(t_disp[2:75], Rᵤᵤ[2:75,100]./Rᵤᵤ[1,100],
    label="Generated Field", xlabel="UΔtM⁻¹", xlim=(0.1,200), ylim=(0.001,10),
    ylabel="R₁₁ (x₀M⁻¹ = 42, 0, 0; ΔxM⁻¹)",legendfontsize=6, xscale=:log10,
    yscale=:log10)
scatter!(Ruu_plot1, t_data, R11_508[:,2],label="Comte-Bellot and Corrsin (1971)")
display(Ruu_plot1)
savefig(Ruu_plot1, string("Ruu_plot1.pdf"))

# T_plot = heatmap(x1, x2, T_xy[:,:,2])
# display(T_plot)
#
T_plot = heatmap(x, y, T_xy[:,:,2], xlabel="x [cm]", ylabel="y [cm]",
    colorbar_title = "T [K]")
display(T_plot)
savefig(T_plot, "T_plot.pdf")

ux_plot = heatmap(x, y, ux_xy[:,:,2], xlabel="x [cm]", ylabel="y [cm]",
    colorbar_title = "uₓ [cm/s]")
display(ux_plot)
savefig(ux_plot, "ux_plot.pdf")


# Derivative of F_11 to obtain E
avg_ux = mean(ux_xy)
var_ux = var(ux_xy)

L_int = trapz(x[1:75],Rᵤᵤ[1:75])
