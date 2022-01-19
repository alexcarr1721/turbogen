using DelimitedFiles
using FFTW
using Interpolations
using HDF5
using SpecialFunctions
using Distributions
using Statistics
using LinearAlgebra
using LaTeXStrings

# Functions
function autocorrelation(
    f::AbstractArray{<:Number,N} where N, # Signal
    x::AbstractArray{<:Number,N} where N,  # Domain
    σ::Number                             # Standard deviation
    )

    # x must be evenly spaced
    dx = x[2] - x[1]

    # Compute mean and variance from Statistics
    # f_mean = Statistics.mean(f)
    # f_var  = Statistics.var(f)
    f_var = σ^2

    # Compute autocorrelation

    # zero pad by length of array
    n = size(f,1)
    fnew = zeros(typeof(f[1]), n*2)
    fnew[1:n] = f #.- f_mean

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
    dfdx = zeros(Float32, size(f,1))
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
        if ( i == size(f,1)-1 )
            index_zero = i + 1
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

# function fg_model(
#     r::AbstractArray{<:Number,N} where N,
#     σ::Number,  # Standard deviation
#     Lo::Number  # Outer length scale
#     )
#
#     f_model = zero(r)
#     g_model = zero(r)
#
#     for i ∈ 1:size(r,1)
#         f_model[i] = σ^2 * (0.5925485*((r[i]/Lo)^(1/3))*besselk(1/3,r[i]/Lo))
#         g_model[i] = σ^2 * 0.5925485 * ((r[i]/Lo)^(1/3))*(
#             (4/3)*besselk(1/3,r[i]/Lo) - 0.5*(r[i]/Lo)*besselk(4/3,r[i]/Lo) )
#     end
#
#     return f_model, g_model
# end

function fg_model(
    r::AbstractArray{<:Number,N} where N,
    σ::Number,  # Standard deviation
    Lo::Number  # Outer length scale
    )

    f_model = zero(r)
    g_model = zero(r)

    for i = 1:size(r,1)
        if r[i] == 0
            f_model[i] = 1.0
            g_model[i] = 1.0
        else
            f_model[i] = 1.0 * (0.5925485*((r[i]/Lo)^(1/3))*besselk(1/3,r[i]/Lo))
            g_model[i] = 1.0 * 0.5925485 * ((r[i]/Lo)^(1/3))*(
                (4/3)*besselk(1/3,r[i]/Lo) - 0.5*(r[i]/Lo)*besselk(4/3,r[i]/Lo) )
        end
    end

    return f_model, g_model
end


# Parameters ###################################################################
L_f = 0.0240 # Longitudinal integral scale
cbc_std = 0.222  # standard deviation of velocity
u0 = 10.0    # 10 meter/sec freestream velocity
################################################################################

# Load data ####################################################################
filename = "cbc.h5"
x1  = h5read(filename, "grid/x1")
x2  = h5read(filename, "grid/x2")
u1  = h5read(filename, "turbulence/u1")
u2  = h5read(filename, "turbulence/u2")
temperature  = h5read(filename, "turbulence/temperature")
q  = h5read(filename, "turbulence/q")
################################################################################

# Scalar Statistics ############################################################
u1_mean = mean(u1)
u2_mean = mean(u2)
u1_std = std(u1)
u2_std = std(u2)
u1_skw = moment(u1, 3)
u2_skw = moment(u2, 3)
u1_kur = moment(u1, 4)
u2_kur = moment(u2, 4)

# Velocity trace <uᵢuᵢ>
vtrace = 0.5*(u1_std^2 + u2_std^2)
################################################################################

# Correlations #################################################################
u1_autocorrelations = zeros(Float32, size(u1,1), size(u1,2), size(u1,3))
u2_autocorrelations = zeros(Float32, size(u1,1), size(u1,2), size(u1,3))

for k = 1:size(u1,3)
    for j = 1:size(u1,2)
        u1_autocorrelations[:,j,k] = autocorrelation(u1[:,j,k], x1, u1_std)
        u2_autocorrelations[:,j,k] = autocorrelation(u2[:,j,k], x1, u2_std)
    end
end

Longitudinal_correlation = zeros(Float32, size(u1,1))
Lateral_correlation = zeros(Float32, size(u1,1))

for i = 1:size(u1,1)
    Longitudinal_correlation[i] = mean(u1_autocorrelations[i,:,:])
    Lateral_correlation[i] = mean(u2_autocorrelations[i,:,:])
end

f_model, g_model = fg_model(x1.-x1[1], cbc_std, L_f/0.7468)
################################################################################

# Taylor microscale ############################################################
fpp = d2fdx2_6(Longitudinal_correlation, x1)
gpp = d2fdx2_6(Lateral_correlation, x1)

λ = sqrt(-1.0/fpp[1])
mean_square_vorticity = 15*vtrace/(λ^2)

Lp = integrate_to_zero(Longitudinal_correlation, x1)
Ln = integrate_to_zero(Lateral_correlation, x1)

# # Write to file
# io = open("turbulent_field.txt", "w")
# write(io, "uf^2     ug^2    Lf    Lg    Taylor Microscale   Mean-Square Vorticity")
# write(io,"\n")
# write(io, string(ux1_std^2, "   ",ux2_std^2,"   ",Lp,"   ",
#     Ln,"   ",λ,"   ",mean_square_vorticity))
# write(io,"\n")
# write(io,"\n")
# write(io,"\n")
# write(io, "Mean ux1: "); write(io, string(ux1_mean))
# write(io,"\n")
# write(io, "Mean ux2: "); write(io, string(ux2_mean))
# write(io,"\n")
# write(io, "Mean ux3: "); write(io, string(ux3_mean))
# write(io,"\n")
# write(io, "Std dev ux1: "); write(io, string(ux1_std))
# write(io,"\n")
# write(io, "Std dev ux2: "); write(io, string(ux2_std))
# write(io,"\n")
# write(io, "Std dev ux3: "); write(io, string(ux3_std))
# write(io,"\n")
# write(io, "1/3(ux1_var + ux2_var + ux3_var): "); write(io, string(vtrace))
# close(io)
################################################################################

# Write to file ################################################################
fid = h5open("cbc_statistics.h5", "w")
fid["grid/x1"] = x1
fid["grid/x2"] = x2
fid["scalar/u1_mean"] = u1_mean
fid["scalar/u2_mean"] = u2_mean
fid["scalar/u1_std"] = u1_std
fid["scalar/u2_std"] = u2_std
fid["scalar/turbulence_intesity"] = vtrace
fid["scalar/Taylor_microscale"] = λ
fid["scalar/mean_square_vorticity"] = mean_square_vorticity
fid["scalar/Lf"] = Lp
fid["scalar/Lg"] = Ln
fid["correlations/f_model"] = f_model
fid["correlations/g_model"] = g_model
fid["correlations/Longitudinal_correlation"] = Longitudinal_correlation
fid["correlations/Lateral_correlation"] = Lateral_correlation
fid["fields/u1"] = u1[:,:,1]
fid["fields/u2"] = u2[:,:,1]
fid["cbc/L_f"] = L_f
fid["cbc/std"] = cbc_std
fid["cbc/u0"] = u0
close(fid)
################################################################################
