#!/usr/bin/env python3

from mpi4py import MPI
from numpy.core.function_base import linspace
from numpy.core.numeric import full
from scipy import special
from scipy.linalg import eig, eigh
from scipy.sparse.linalg import eigs, eigsh
import numpy
import h5py
import sys
import time

def most(z, ustar, Tstar, qstar, Lobukhov, z0, zr, Tr, qr):
    kappa = 0.4
    Pt = 0.95
    Gammad = 0.0098

    n = len(z)
    v = numpy.zeros((n))
    T = numpy.zeros((n))
    q = numpy.zeros((n))

    for i in range(0,n): # Although it might seem like this will loop
        # from 0 to 100, it will actually go from 0 to 99. Python is wierd.
        if z[i] == 0:
            ztemp = z0/1.5
            v[i] = (ustar/kappa)*(numpy.log(ztemp/z0) \
                - wilson_psi(ztemp/Lobukhov, 2) + wilson_psi(z0/Lobukhov, 2))
            T[i] = Tr - (ztemp - zr)*Gammad + (Pt*Tstar/kappa)*( numpy.log(ztemp/zr) \
                - wilson_psi(ztemp/Lobukhov, 1) + wilson_psi(z0/Lobukhov, 1) )
            q[i] = qr - (ztemp - zr)*Gammad + (Pt*qstar/kappa)*( numpy.log(ztemp/zr) \
                - wilson_psi(ztemp/Lobukhov, 1) + wilson_psi(z0/Lobukhov, 1) )
        else:
            v[i] = (ustar/kappa)*(numpy.log(z[i]/z0) \
                - wilson_psi(z[i]/Lobukhov, 2) + wilson_psi(z0/Lobukhov, 2))
            T[i] = Tr - (z[i] - zr)*Gammad + (Pt*Tstar/kappa)*( numpy.log(z[i]/zr) \
                - wilson_psi(z[i]/Lobukhov, 1) + wilson_psi(z0/Lobukhov, 1) )
            q[i] = qr - (z[i] - zr)*Gammad + (Pt*qstar/kappa)*( numpy.log(z[i]/zr) \
                - wilson_psi(z[i]/Lobukhov, 1) + wilson_psi(z0/Lobukhov, 1) )
        
    
    return v, T, q

def wilson_psi(zeta, m):
    if m == 1: # Temperature
        if zeta < 0:
            psi = 2.0*numpy.log( (1 + numpy.sqrt(1 + \
                7.9*(numpy.abs(zeta)**(2/3)) ) )/2 )
        else:
            psi = -8.4*zeta 
    elif m == 2: # velocity
        if zeta < 0:
            psi = 2.0*numpy.log( (1 + numpy.sqrt(1 + \
                3.6*(numpy.abs(zeta)**(2/3)) ) )/2 )
        else:
            psi = -5.3*zeta 
    else: # m = 1 by default
        if zeta < 0:
            psi = 2.0*numpy.log( (1 + numpy.sqrt(1 + \
                7.9*(numpy.abs(zeta)**(2/3)) ) )/2 )
        else:
            psi = -8.4*zeta 

    return psi 

def temperature_statistics(z, Lobukhov, Tstar):
    n = len(z)
    sigmaT = numpy.zeros((n))
    LT = numpy.zeros((n))
    min_index = numpy.abs(z - 200.0).argmin()

    for i in range(0,n):
        if z[i] > z[min_index]:
            sigmaT[i] = sigmaT[min_index]
            LT[i] = LT[min_index]
        else: 
            sigmaT[i] = numpy.sqrt( 4.0*(Tstar**2)/( \
                (1 + 10*(numpy.abs(z[i]/Lobukhov)))**(2/3) ) )
            LT[i] = 2.0*z[i]*(1 + 7*(numpy.abs(z[i]/Lobukhov)))/(1 \
                + 10*(numpy.abs(z[i]/Lobukhov)))

    return sigmaT, LT

def shear_statistics(z, Lobukhov, ustar):
    n = len(z)
    sigmaS = numpy.zeros((n))
    LS = numpy.zeros((n))
    min_index = numpy.abs(z - 200.0).argmin()

    for i in range(0,n):
        if z[i] > z[min_index]:
            sigmaS[i] = sigmaS[min_index]
            LS[i] = LS[min_index]
        else: 
            sigmaS[i] = numpy.sqrt(3.0*(ustar**2))
            LS[i] = 1.8*z[i]

    return sigmaS, LS

def bouyancy_statistics(z, Lobukhov, wstar):
    n = len(z)
    sigmaB = numpy.zeros((n))
    LB = numpy.zeros((n))

    for i in range(0,n):
        sigmaB[i] = numpy.sqrt(0.35*(wstar**2))
        LB[i] = 0.23*z[i]

    return sigmaB, LB

def temperature_correlation(z, z0, LT, sigmaT, kx):
    if LT == 0:
        zeta = 5000.0
    else:
        zeta = (numpy.abs(z - z0)/LT)*numpy.sqrt(1.0 + (kx**2 * LT**2))


    if zeta == 0:
        BT = (2*(sigmaT**2)*LT/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            (1/(1 + (kx**2 * LT**2)))**(5/6))*( \
            (11/3)*special.gamma(5/6)/2)
    elif zeta > 750.0:
        BT = 0.0
    else:
        BT = (2*(sigmaT**2)*LT/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            ((zeta/2)/(1 + (kx**2 * LT**2)))**(5/6))*( \
            (11/3)*special.kv(5/6,zeta))

    return BT

def lateral_correlation(z, z0, LV, sigmaV, kx):
    if LV == 0:
        zeta = 5000.0
    else:
        zeta = (numpy.abs(z - z0)/LV)*numpy.sqrt(1.0 + (kx**2 * LV**2))

    if zeta == 0:
        B22 = (2*(sigmaV**2)*LV/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            (1/(1 + (kx**2 * LV**2)))**(5/6))*( \
            (4/3)*special.gamma(5/6)/2 - \
            (1/(1 + (kx**2 * LV**2)))*special.gamma(11/6)/2)
    elif zeta > 1500.0:
        B22 = 0.0
    else:
        B22 = (2*(sigmaV**2)*LV/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            ((zeta/2)/(1 + (kx**2 * LV**2)))**(5/6))*( \
            (4/3)*special.kv(5/6,zeta) - \
            ((zeta/2)/(1 + (kx**2 * LV**2)))*special.kv(11/6,zeta) + \
            (zeta/2)*special.kv(1/6,zeta))

    return B22

def longitudinal_correlation(z, z0, LV, sigmaV, kx):
    if LV == 0:
        zeta = 5000.0
    else:
        zeta = (numpy.abs(z - z0)/LV)*numpy.sqrt(1.0 + (kx**2 * LV**2))

    if zeta == 0:
        B11 = (2*(sigmaV**2)*LV/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            (1/(1 + (kx**2 * LV**2)))**(5/6))*( \
            special.gamma(5/6)/2)
    elif zeta > 1500.0:
        B11 = 0.0
    else:
        B11 = (2*(sigmaV**2)*LV/(numpy.sqrt(numpy.pi)*special.gamma(1/3)))*( \
            ((zeta/2)/(1 + (kx**2 * LV**2)))**(5/6))*( \
            special.kv(5/6,zeta) - \
            (zeta/2)*special.kv(1/6,zeta))

    return B11

def frequency_bins(x):
    f = numpy.zeros_like(x)
    N = len(x)
    fs = N / (abs(x[N-1] - x[0]))
    if N%2 == 0:
        f = numpy.arange(-N/2,N/2) * (fs / N)
    else:
        f = numpy.arange(-(N - 1)/2,(N)/2) * (fs / N)

    f = numpy.fft.fftshift(f)
    return f

def power_iteration(A, Omega, power_iter = 3):
    Y = A @ Omega
    for q in range(power_iter):
        Y = A @ (A.T @ Y)
    Q, _ = numpy.linalg.qr(Y)
    return Q 
    
def rsvd(A, Omega):
    Q = power_iteration(A, Omega)
    B = Q.T @ A
    u_tilde, s, v = numpy.linalg.svd(B, full_matrices = 0)
    u = Q @ u_tilde
    return u, s, v

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
proc = comm.Get_rank()

# Parameters
kappa = 0.4
Pt = 0.95
Gammad = 0.0098
gravity = 9.81
Tr = 293.15     # Temperature at zr
qr = 0.0        # Humidity at observation height
zr = 10.0       # Observation height
z0 = 0.01       # Roughness length
zi = 1000.0     # ABL height
c0 = 343.0      # Speed of sound
Q0 = 0.01       # Sensible surface heat flux
ustar = 0.05    # Shear velocity
Lobukhov = -1.0*(ustar**3)*Tr/(kappa*gravity*Q0)
Tstar = -Q0/ustar
wstar = (zi*gravity*Q0/Tr)**(1/3)    
qstar = 0.0     # Humidity scale (dry air for now)

# Domain #################################################
N = 64
M = 64
modes = 50
local_size = int(N/nproc)
x = numpy.linspace(0, 200.0, N)
z = numpy.linspace(0, 200.0, M)
fx = frequency_bins(x)
kx = 2*numpy.pi*fx
##########################################################

# Allocate matrices ######################################
RTT = numpy.zeros((M, M), dtype=numpy.float32)
R11 = numpy.zeros((M, M), dtype=numpy.float32)
R22 = numpy.zeros((M, M), dtype=numpy.float32)
UTT = numpy.zeros((M, modes, local_size), dtype=numpy.float32)
U11 = numpy.zeros((M, modes, local_size), dtype=numpy.float32)
U22 = numpy.zeros((M, modes, local_size), dtype=numpy.float32)
VTT = numpy.zeros((modes, M, local_size), dtype=numpy.float32)
V11 = numpy.zeros((modes, M, local_size), dtype=numpy.float32)
V22 = numpy.zeros((modes, M, local_size), dtype=numpy.float32)
STT = numpy.zeros((modes, local_size), dtype=numpy.float32)
S11 = numpy.zeros((modes, local_size), dtype=numpy.float32)
S22 = numpy.zeros((modes, local_size), dtype=numpy.float32)

# Iterate to find real Lobukhov and wstar
for i in range(1,4):
    v, T, q = most(z, ustar, Tstar, qstar, Lobukhov, z0, zr, Tr, qr)
    Lobukhov = -1.0*(ustar**3)*T[0]/(kappa*gravity*Q0)
    wstar = (zi*gravity*Q0/T[0])**(1/3)

# Find statistics in the vertical direction
sigmaT, LT = temperature_statistics(z, Lobukhov, Tstar)
sigmaS, LS = shear_statistics(z, Lobukhov, ustar)
sigmaB, LB = bouyancy_statistics(z, Lobukhov, wstar)

# Find all of the eigenvalues and eigenvectors ########################
for k in range(0,local_size):
    global_ind = k + proc*local_size
    for j in range(0,M):
        for i in range(0,M):
            RTT[i,j] = temperature_correlation(z[i], z[j], LT[j], sigmaT[j], kx[global_ind])
            R11[i,j] = longitudinal_correlation(z[i], z[j], LS[j], sigmaS[j], kx[global_ind]) \
                + longitudinal_correlation(z[i], z[j], LB[j], sigmaB[j], kx[global_ind])
            R22[i,j] = lateral_correlation(z[i], z[j], LS[j], sigmaS[j], kx[global_ind]) \
                + lateral_correlation(z[i], z[j], LB[j], sigmaB[j], kx[global_ind])

    omegaT = numpy.random.randn(RTT.shape[1], modes)
    UTT[:,:,k], STT[:,k], VTT[:,:,k] = rsvd(RTT, omegaT)
    omega1 = numpy.random.randn(R11.shape[1], modes)
    U11[:,:,k], S11[:,k], V11[:,:,k] = rsvd(R11, omega1)
    omega2 = numpy.random.randn(R22.shape[1], modes)
    U22[:,:,k], S22[:,k], V22[:,:,k] = rsvd(R22, omega2)

########################################################################

# for j in range(0,M):
#     for i in range(0,M):
#         RTT[i,j] = temperature_correlation(z[i], z[j], LT[j], sigmaT[j], kx[0])
#         R11[i,j] = longitudinal_correlation(z[i], z[j], LS[j], sigmaS[j], kx[0]) \
#             + longitudinal_correlation(z[i], z[j], LB[j], sigmaB[j], kx[0])
#         R22[i,j] = lateral_correlation(z[i], z[j], LS[j], sigmaS[j], kx[0]) \
#             + lateral_correlation(z[i], z[j], LB[j], sigmaB[j], kx[0])

# start = time.time()
# u, s, v = numpy.linalg.svd(RTT, full_matrices= 0 )
# end = time.time()
# print ("Time elapsed for full svd:", end - start)


# start = time.time()
# omega = numpy.random.randn(RTT.shape[1], modes)
# unew, snew, vnew = rsvd(RTT, omega)
# end = time.time()
# print ("Time elapsed for randomized svd:", end - start)

# L2norm = numpy.sqrt(numpy.sum(numpy.abs(snew - s[:10])**2))
# print(L2norm)

# Write to file ########################################################
f = h5py.File('svd.h5', 'w', driver='mpio', comm=MPI.COMM_WORLD)

dset = f.create_dataset('temperature/sigma', (modes, M), dtype=numpy.float32)
dset[:,(0 + proc*local_size):((proc+1)*local_size)] = STT
dset = f.create_dataset('u1/sigma', (modes, M), dtype=numpy.float32)
dset[:,(0 + proc*local_size):((proc+1)*local_size)] = S11
dset = f.create_dataset('u2/sigma', (modes, M), dtype=numpy.float32)
dset[:,(0 + proc*local_size):((proc+1)*local_size)] = S22
dset = f.create_dataset('temperature/U', (M, modes, M), dtype=numpy.float32)
dset[:,:,(0 + proc*local_size):((proc+1)*local_size)] = UTT
dset = f.create_dataset('u1/U', (M, modes, M), dtype=numpy.float32)
dset[:,:,(0 + proc*local_size):((proc+1)*local_size)] = U11
dset = f.create_dataset('u2/U', (M, modes, M), dtype=numpy.float32)
dset[:,:,(0 + proc*local_size):((proc+1)*local_size)] = U22

dset = f.create_dataset('stddev/temperature', (M), dtype=numpy.float32)
dset = sigmaT
dset = f.create_dataset('stddev/shear', (M), dtype=numpy.float32)
dset = sigmaS
dset = f.create_dataset('stddev/bouyancy', (M), dtype=numpy.float32)
dset = sigmaB
dset = f.create_dataset('vonkarmanscale/temperature', (M), dtype=numpy.float32)
dset = LT
dset = f.create_dataset('vonkarmanscale/shear', (M), dtype=numpy.float32)
dset = LS
dset = f.create_dataset('vonkarmanscale/bouyancy', (M), dtype=numpy.float32)
dset = LB

f.close()
########################################################################