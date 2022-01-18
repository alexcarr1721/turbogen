################################################################################
#                                                                              #
#  Generate Inhomogeneous turbulence                                           #
#                                                                              #
#   Author:     Alexander Carr                                                 #
#   Title:      Ph.D. Candidate, University of Florida                         #
#   Date:       05/26/2021                                                     #
#   Contact:    alexcarr.1721@gmail.com                                        #
#                                                                              #
################################################################################

# Load Libraries ###############################################################
using LinearAlgebra
using Interpolations
using Statistics
using Distributions
using SpecialFunctions
# using LowRankApprox
using RandomizedLinAlg
using FFTW
using HDF5
using Plots
gr()
################################################################################

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
    elseif ( ζ > 2000 )
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

# Set up domain ################################################################
N = 4096
M = 4096
modes = 50
x₁ = range(0,stop=4000.0,length=N)
x₃ = range(0,stop=1500.0,length=M)
f₁ = fourier_space(x₁)
k₁ = 2*π*f₁
################################################################################

# Allocate cross spectrum, ritz values, ritz vectors ###########################
sT  = zeros(ComplexF32, N, M)
s1  = zeros(ComplexF32, N, M)
s3  = zeros(ComplexF32, N, M)
T  = zeros(Float32, N, M)
u1  = zeros(Float32, N, M)
u3  = zeros(Float32, N, M)
U₀₀ = zeros(Float32, M, modes)
U₁₁ = zeros(Float32, M, modes)
U₃₃ = zeros(Float32, M, modes)
S₀₀ = zeros(Float32, modes)
S₁₁ = zeros(Float32, modes)
S₃₃ = zeros(Float32, modes)
################################################################################

# Make directory ###############################################################
try
    mkdir("turbulence")
catch
    # do nothing
end
################################################################################

# Generate turbulent field #####################################################
read_filename = "svd.h5"
fid = h5open(read_filename, "r")
dset_stdtemperature = fid["stddev/temperature"]
σT = dset_stdtemperature[:]
dset_stdshear = fid["stddev/shear"]
σs = dset_stdshear[:]
dset_stdbouyancy = fid["stddev/bouyancy"]
σb = dset_stdbouyancy[:]
dset_vonkarmanlengthtemperature = fid["vonkarmanlength/temperature"]
LT = dset_vonkarmanlengthtemperature[:]
dset_vonkarmanlengthshear = fid["vonkarmanlength/shear"]
Ls = dset_vonkarmanlengthshear[:]
dset_vonkarmanlengthbouyancy = fid["vonkarmanlength/bouyancy"]
Lb = dset_vonkarmanlengthbouyancy[:]
for j = 1:N
    # Read singular values and vectors from file ###############################
    dset_u1U = fid["u1/U"]
    U₁₁ = dset_u1U[:,:,j]
    dset_u2U = fid["u2/U"]
    U₃₃ = dset_u2U[:,:,j]
    dset_TU = fid["temperature/U"]
    U₀₀ = dset_TU[:,:,j]
    dset_u1sigma = fid["u1/sigma"]
    S₁₁ = dset_u1sigma[:,j]
    dset_u2sigma = fid["u2/sigma"]
    S₃₃ = dset_u2sigma[:,j]
    dset_Tsigma = fid["temperature/sigma"]
    S₀₀ = dset_Tsigma[:,j]
    ############################################################################

    # Compute array of random numbers
    γ₀₀ = rand(Uniform(-1.73205,1.73205), modes)*(1 + 1im)
    γ₁₁ = rand(Uniform(-1.73205,1.73205), modes)*(1 + 1im)
    γ₃₃ = rand(Uniform(-1.73205,1.73205), modes)*(1 + 1im)

    # Summation
    for i ∈ 1:M
        sT[j,i] = 1.0*sum(sqrt.(S₀₀).*γ₀₀.*U₀₀[i,:])
        s1[j,i] = 1.0*sum(sqrt.(S₁₁).*γ₁₁.*U₁₁[i,:])
        s3[j,i] = 1.0*sum(sqrt.(S₃₃).*γ₃₃.*U₃₃[i,:])
    end
end
close(fid)

# Inverse FFT
T = real(ifft(sT,1))*N
u1 = real(ifft(s1,1))*N
u3 = real(ifft(s3,1))*N

# Compute variance at each z location and re-adjust to actual variance
u1std = std(u1)
u3std = std(u3)
Tstd = zeros(Float64, M)
for j = 1:M
    for i = 1:N
        u1[i,j] = ((σs[j]+σb[j])/u1std)*u1[i,j]
        u3[i,j] = ((σs[j]+σb[j])/u3std)*u3[i,j]
    end
    Tstd[j] = std(T[:,j])
end

# Compare average Tstd to average σT
for j = 1:M
    for i = 1:N
        T[i,j] = (mean(σT)/mean(Tstd))*T[i,j]
    end
end

# Write to hdf5 file
directory = "turbulence"
filename = "field"
run_number = ARGS[1]
fid = h5open(string(directory,"/",filename,run_number,".h5"), "w")
fid["turbulence/u1"] = u1
fid["turbulence/u2"] = u3
fid["turbulence/temperature"] = T
fid["stddev/temperature"] = σT
fid["stddev/shear"] = σs
fid["stddev/bouyancy"] = σb
fid["vonkarmanlength/temperature"] = LT
fid["vonkarmanlength/shear"] = Ls
fid["vonkarmanlength/bouyancy"] = Lb
close(fid)
