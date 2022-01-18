################################################################################
#                                                                              #
#     Gather data from output files and save to a single hdf5 file             #
#                                                                              #
################################################################################

using DelimitedFiles
using Interpolations
using HDF5

count = 2 # number of datasets

# Step 1) Read grid and other parameters that are constant for each dataset ####
filename_starter = "svd/svd_k"
sigmaT = h5read(string(filename_starter,1,".h5"),"temperature/sigma")
sigmau1 = h5read(string(filename_starter,1,".h5"),"u1/sigma")
sigmau2 = h5read(string(filename_starter,1,".h5"),"u2/sigma")
UT = h5read(string(filename_starter,1,".h5"),"temperature/U")
Uu1 = h5read(string(filename_starter,1,".h5"),"u1/U")
Uu2 = h5read(string(filename_starter,1,".h5"),"u2/U")
stddev_temperature = h5read(string(filename_starter,1,".h5"),"stddev/temperature")
stddev_shear = h5read(string(filename_starter,1,".h5"),"stddev/shear")
stddev_bouyancy = h5read(string(filename_starter,1,".h5"),"stddev/bouyancy")
vonkarmanlength_temperature = h5read(string(filename_starter,1,".h5"),"vonkarmanlength/temperature")
vonkarmanlength_shear = h5read(string(filename_starter,1,".h5"),"vonkarmanlength/shear")
vonkarmanlength_bouyancy = h5read(string(filename_starter,1,".h5"),"vonkarmanlength/bouyancy")
################################################################################

filename = "svd.h5"
fid = h5open(filename, "w")
fid["stddev/temperature"] = stddev_temperature
fid["stddev/shear"] = stddev_shear
fid["stddev/bouyancy"] = stddev_bouyancy
fid["vonkarmanlength/temperature"] = vonkarmanlength_temperature
fid["vonkarmanlength/shear"] = vonkarmanlength_shear
fid["vonkarmanlength/bouyancy"] = vonkarmanlength_bouyancy
n = size(UT,1)
modes = size(UT,2)
dset_sigmaT = create_dataset(fid, "temperature/sigma", datatype(Float32),
    dataspace(modes,count), chunk=(modes,1), compress=3)
dset_sigma1 = create_dataset(fid, "u1/sigma", datatype(Float32),
    dataspace(modes,count), chunk=(modes,1), compress=3)
dset_sigma2 = create_dataset(fid, "u2/sigma", datatype(Float32),
    dataspace(modes,count), chunk=(modes,1), compress=3)
dset_UT = create_dataset(fid, "temperature/U", datatype(Float32),
    dataspace(n,modes,count), chunk=(n,modes,1), compress=3)
dset_U1 = create_dataset(fid, "u1/U", datatype(Float32),
    dataspace(n,modes,count), chunk=(n,modes,1), compress=3)
dset_U2 = create_dataset(fid, "u2/U", datatype(Float32),
    dataspace(n,modes,count), chunk=(n,modes,1), compress=3)
dset_sigmaT[:,1] = sigmaT
dset_sigma1[:,1] = sigmau1
dset_sigma2[:,1] = sigmau2
dset_UT[:,:,1] = UT
dset_U1[:,:,1] = Uu1
dset_U2[:,:,1] = Uu2
if (count > 1)
    for i = 2:count
        sigmaT = h5read(string(filename_starter,i,".h5"),"temperature/sigma")
        sigmau1 = h5read(string(filename_starter,i,".h5"),"u1/sigma")
        sigmau2 = h5read(string(filename_starter,i,".h5"),"u2/sigma")
        UT = h5read(string(filename_starter,i,".h5"),"temperature/U")
        Uu1 = h5read(string(filename_starter,i,".h5"),"u1/U")
        Uu2 = h5read(string(filename_starter,i,".h5"),"u2/U")
        # Save to new file
        dset_sigmaT[:,i] = sigmaT
        dset_sigma1[:,i] = sigmau1
        dset_sigma2[:,i] = sigmau2
        dset_UT[:,:,i] = UT
        dset_U1[:,:,i] = Uu1
        dset_U2[:,:,i] = Uu2
    end
end
close(fid)
