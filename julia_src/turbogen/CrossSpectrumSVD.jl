################################################################################
#                                                                              #
#  Compute the SVD of the von Karman cross spectral model                      #
#                                                                              #
#   Author:     Alexander Carr                                                 #
#   Title:      Ph.D. Candidate, University of Florida                         #
#   Date:       12/14/2021                                                     #
#   Contact:    alexcarr.1721@gmail.com                                        #
#                                                                              #
################################################################################

# Load Libraries ###############################################################
using MPI
using HDF5
using LinearAlgebra
using Interpolations
using Statistics
using Distributions
using SpecialFunctions
# using LowRankApprox
using RandomizedLinAlg
using FFTW
using DelimitedFiles
using Plots
MPI.Init()
gr()
################################################################################

# Functions ####################################################################
function most(
    z::AbstractArray{<:Number,N} where N,
    u_star::Number,
    T_star::Number,
    L_0::Number,
    T_r::Number,
    z_r::Number,
    z_0::Number,
    Γ_d::Number,
    P_t::Number,
    κ::Number)

    T = zeros(Float64, size(z,1))
    v = zeros(Float64, size(z,1))

    for i = 1:size(z,1)
        if z[i] == 0
            ztemp = z_0/1.5
            T[i] = T_r - (ztemp - z_r)*Γ_d + (P_t*T_star/κ)*(log(ztemp/z_r) -
                Wilsonpsi(ztemp/L_0,1) + Wilsonpsi(z_r/L_0,1))
            v[i] = (u_star/κ)*(log(ztemp/z_0) - Wilsonpsi(ztemp/L_0,2) + Wilsonpsi(z_0/L_0,2) )
        else
            T[i] = T_r - (z[i] - z_r)*Γ_d + (P_t*T_star/κ)*(log(z[i]/z_r) -
                Wilsonpsi(z[i]/L_0,1) + Wilsonpsi(z_r/L_0,1))
            v[i] = (u_star/κ)*(log(z[i]/z_0) - Wilsonpsi(z[i]/L_0,2) + Wilsonpsi(z_0/L_0,2) )
        end
    end

    return v, T
end

function Wilsonpsi(ξ::Number, m::Number)
    if m == 1 # h
        if ξ < 0
            Ψ = 2*log( (1 + sqrt(1 + 7.9*(abs(ξ)^(2/3)) ) )/(2) )
        else
            Ψ = -8.4*ξ
        end
    elseif m == 2 # m
        if ξ < 0
            Ψ = 2*log( (1 + sqrt(1 + 3.6*(abs(ξ)^(2/3)) ) )/(2) )
        else
            Ψ = -5.3*ξ
        end
    end

    return Ψ
end

function TemperatureStatistics(
    z::AbstractArray{<:Number,N} where N,
    Lo::Number,
    Tstar::Number
    )

    n = size(z,1)
    σT = zeros(Float32, n)
    LT = zeros(Float32, n)
    min_index = argmin(abs.(z .- 200))

    for i = 1:n
        # Only applicable till about 200 meters off the ground, after that keep constant
        if ( z[i] > z[min_index] )
            σT[i] = σT[min_index]
            LT[i] = LT[min_index]
        else
            σT[i] = sqrt( 4.0*(Tstar^2)/( (1 + 10*(abs.(z[i]/Lo)))^(2/3) ) )
            LT[i] = 2.0*z[i]*(1 + 7*(abs.(z[i]/Lo)))/(1 + 10*(abs.(z[i]/Lo)))
        end
    end

    return σT, LT
end

function ShearStatistics(
    z::AbstractArray{<:Number,N} where N,
    Lo::Number,
    ustar::Number
    )

    n = size(z,1)
    σ = zeros(Float32, n)
    L = zeros(Float32, n)
    min_index = argmin(abs.(z .- 200))

    for i = 1:n
        # Only applicable till about 200 meters off the ground, after that keep constant
        if ( z[i] > z[min_index] )
            σ[i] = σ[min_index]
            L[i] = L[min_index]
        else
            σ[i] = sqrt( 3.0*(ustar^2) )
            L[i] = 1.8*z[i]
        end
    end

    return σ, L
end

function BouyancyStatistics(
    z::AbstractArray{<:Number,N} where N,
    Lo::Number,
    wstar::Number,
    zi::Number
    )

    n = size(z,1)
    σ = zeros(Float32, n)
    L = zeros(Float32, n)

    for i = 1:n
        σ[i] = sqrt( 0.35*(wstar^2) )
        L[i] = 0.23*zi
    end

    return σ, L
end

function TemperatureCorrelation(
    z::Number,
    z0::Number,
    LT::Number,
    σT::Number,
    kx::Number
    )

    ζ = (abs(z - z0)/LT)*sqrt(1 + (kx^2 * LT^2))
    if ( LT == 0 )
        ζ = Inf
    end

    if ( ζ == 0 )
        BT = (2*(σT^2)*LT/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            (1/(1 + (kx^2 * LT^2)))^(5/6))*(
            (11/3)*SpecialFunctions.gamma(5/6)/2)
    elseif ( ζ > 500 )
        BT = 0.0
    else
        BT = (2*(σT^2)*LT/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            ((ζ/2)/(1 + (kx^2 * LT^2)))^(5/6))*(
            (11/3)*SpecialFunctions.besselk(5/6,ζ))
    end

    return BT
end

function LongitudinalCorrelation(
    z::Number,
    z0::Number,
    Lv::Number,
    σv::Number,
    kx::Number
    )

    ζ = (abs(z - z0)/Lv)*sqrt(1 + (kx^2 * Lv^2))
    if ( Lv == 0 )
        ζ = Inf
    end


    if ( ζ == 0 )
        B11 = (2*(σv^2)*Lv/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            (1/(1 + (kx^2 * Lv^2)))^(5/6))*(
            SpecialFunctions.gamma(5/6)/2)
    elseif ( ζ > 2000 )
        B11 = 0.0
    else
        B11 = (2*(σv^2)*Lv/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            ((ζ/2)/(1 + (kx^2 * Lv^2)))^(5/6))*(
            SpecialFunctions.besselk(5/6,ζ) -
            (ζ/2)*SpecialFunctions.besselk(1/6,ζ))
    end

    return B11
end

function LateralCorrelation(
    z::Number,
    z0::Number,
    Lv::Number,
    σv::Number,
    kx::Number
    )

    ζ = (abs(z - z0)/Lv)*sqrt(1 + (kx^2 * Lv^2))
    if ( Lv == 0 )
        ζ = Inf
    end

    if ( ζ == 0 )
        B22 = (2*(σv^2)*Lv/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            (1/(1 + (kx^2 * Lv^2)))^(5/6))*(
            (4/3)*SpecialFunctions.gamma(5/6)/2 -
            (1/(1 + (kx^2 * Lv^2)))*SpecialFunctions.gamma(11/6)/2)
    elseif ( ζ > 2000 )
        B22 = 0.0
    else
        B22 = (2*(σv^2)*Lv/(sqrt(π)*SpecialFunctions.gamma(1/3)))*(
            ((ζ/2)/(1 + (kx^2 * Lv^2)))^(5/6))*(
            (4/3)*SpecialFunctions.besselk(5/6,ζ) -
            ((ζ/2)/(1 + (kx^2 * Lv^2)))*SpecialFunctions.besselk(11/6,ζ) +
            (ζ/2)*SpecialFunctions.besselk(1/6,ζ))
    end

    return B22
end

function fourier_space(x::AbstractArray{<:Number,N} where {N})

    f = zero(x)
    N = size(x, 1)
    fs = N / (abs(x[N] - x[1]))
    if mod(N, 2) == 0
        f = (-N/2:N/2-1) * (fs / N)
    else
        f = (-(N - 1)/2:(N-1)/2) * (fs / N)
    end

    f = fftshift(f)
    return f
end
################################################################################

# Communication variables ######################################################
comm = MPI.COMM_WORLD
info = MPI.Info()
proc = MPI.Comm_rank(comm)
nproc = MPI.Comm_size(comm)
################################################################################

# Read input file and weather file #############################################
# input_file_read = readdlm("input.inp")
# weather_file_read = readdlm("weather.inp")
################################################################################

# Set up domain ################################################################
N = 4096 #input_file_read[28,1]
M = 4096 #input_file_read[29,1]
modes = 50
local_size = Int(N/nproc)
x₁ = range(0.0,stop=4000.0,length=N)
x₃ = range(0.0,stop=1500.0,length=M)
f₁ = fourier_space(x₁)
k₁ = 2*π*f₁
################################################################################

# Allocate cross spectrum, ritz values, ritz vectors ###########################
Θ₀₀ = zeros(Float64, M, M)
Θ₁₁ = zeros(Float64, M, M)
Θ₃₃ = zeros(Float64, M, M)
U₀₀ = zeros(Float64, M, modes, local_size)
U₁₁ = zeros(Float64, M, modes, local_size)
U₃₃ = zeros(Float64, M, modes, local_size)
S₀₀ = zeros(Float64, modes, local_size)
S₁₁ = zeros(Float64, modes, local_size)
S₃₃ = zeros(Float64, modes, local_size)
################################################################################

# Weather parameters for the turbulence ########################################
Tr = 293.15 # weather_file_read[11,1]
κ = 0.4
Pt = 0.95
Gamma_d = 0.0098
zr = 10.0 #weather_file_read[12,1]
z0 = 0.01 #weather_file_read[13,1]
zi = 1000.0 #weather_file_read[14,1]
c0 = 343.0
Q0 = 0.01 #weather_file_read[10,1]
ustar = 0.05 #weather_file_read[6,1]
Lobukhov = -1.0*(ustar^3)*Tr/(κ*9.81*Q0) #weather_file_read[9,1] #-1.0*(ustar^3)*Tr/(κ*9.81*Q0)
Tstar = -Q0/ustar #weather_file_read[8,1] #-Q0/ustar
wstar = (zi*9.81*Q0/Tr)^(1/3) #weather_file_read[7,1] #(zi*9.81*Q0/Tr)^(1/3)
for i = 1:3
    global Lobukhov, Tstar, ustar, Q0, zr, z0, Tr, Gamma_d, Pt, κ, wstar
    v, T = most(x₃, ustar,Tstar, Lobukhov, Tr, zr, z0, Gamma_d, Pt, κ)
    # Update Lobukhov
    Lobukhov = -1.0*(ustar^3)*T[1]/(κ*9.81*Q0)
    wstar = (zi*9.81*Q0/T[1])^(1/3)
end

σT, LT = TemperatureStatistics(x₃, Lobukhov, Tstar)
σs, Ls = ShearStatistics(x₃, Lobukhov, ustar)
σb, Lb = BouyancyStatistics(x₃, Lobukhov, wstar, zi)
################################################################################

# Compute the singular values and vectors for each cross-spectrum ##############
for k = 1:local_size
    global_ind = k + proc*local_size
    # Construct spectrum/correlation function ##################################
    for j = 1:M
        for i = 1:M
            Θ₀₀[i,j] = TemperatureCorrelation(x₃[i],x₃[j],LT[j],σT[j],k₁[global_ind])
            Θ₁₁[i,j] = LongitudinalCorrelation(x₃[i],x₃[j],Ls[j],σs[j],k₁[global_ind]) +
                LongitudinalCorrelation(x₃[i],x₃[j],Lb[j],σb[j],k₁[global_ind])          # Incorporates both shear and bouyancy
            Θ₃₃[i,j] = LateralCorrelation(x₃[i],x₃[j],Ls[j],σs[j],k₁[global_ind]) +
                LateralCorrelation(x₃[i],x₃[j],Lb[j],σb[j],k₁[global_ind])
        end
    end
    ############################################################################

    # Perfom SVD ###############################################################
    F₀₀ = rsvd_fnkz(Θ₀₀, 50)
    F₁₁ = rsvd_fnkz(Θ₁₁, 50)
    F₃₃ = rsvd_fnkz(Θ₃₃, 50)
    ############################################################################

    # Save to S and U ##########################################################
    U₀₀[:,:,k] = F₀₀.U
    S₀₀[:,k] = F₀₀.S
    U₁₁[:,:,k] = F₁₁.U
    S₁₁[:,k] = F₁₁.S
    U₃₃[:,:,k] = F₃₃.U
    S₃₃[:,k] = F₃₃.S
    ############################################################################
end
################################################################################

# Write in parallel to HDF5 file ###############################################
filename = "svd.h5" #
fid = h5open(filename, "w", comm, info)
# Temperature cross spectrum modes
global_dims = (modes, N)
dset = create_dataset(fid, "temperature/sigma", datatype(eltype(S₀₀)), dataspace(global_dims), chunk=(modes, N/nproc))
dset[:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = S₀₀
# Temperature cross spectrum singular vector
global_dims = (M, modes, N)
dset = create_dataset(fid, "temperature/U", datatype(eltype(U₀₀)), dataspace(global_dims), chunk=(M, modes, N/nproc))
dset[:,:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = U₀₀
# Longitudinal velocity cross spectrum modes
global_dims = (modes, N)
dset = create_dataset(fid, "u1/sigma", datatype(eltype(S₁₁)), dataspace(global_dims), chunk=(modes, N/nproc))
dset[:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = S₁₁
# Longitudinal velocity cross spectrum singular vector
global_dims = (M, modes, N)
dset = create_dataset(fid, "u1/U", datatype(eltype(U₁₁)), dataspace(global_dims), chunk=(M, modes, N/nproc))
dset[:,:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = U₁₁
# Vertical velocity cross spectrum modes
global_dims = (modes, N)
dset = create_dataset(fid, "u2/sigma", datatype(eltype(S₃₃)), dataspace(global_dims), chunk=(modes, N/nproc))
dset[:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = S₃₃
# Vertical velocity cross spectrum singular vector
global_dims = (M, modes, N)
dset = create_dataset(fid, "u2/U", datatype(eltype(U₃₃)), dataspace(global_dims), chunk=(M, modes, N/nproc))
dset[:,:,Int(1+proc*(N/nproc)):Int((proc+1)*(N/nproc))] = U₃₃
# Case number and Weather data
fid["stddev/temperature"] = σT
fid["stddev/shear"] = σs
fid["stddev/bouyancy"] = σb
fid["integralscale/temperature"] = LT
fid["integralscale/shear"] = Ls
fid["integralscale/bouyancy"] = Lb
close(fid)
################################################################################
