using HDF5

filename = "cbc.h5"
num_of_fields = 20

x1 = h5read(string("field",1,".h5"),"x1")
x2 = h5read(string("field",1,".h5"),"x2")

fid = h5open(filename, "w")
n = size(x1,1)
m = size(x2,1)
fid["grid/x1"] = x1
fid["grid/x2"] = x2
dset_u1 = create_dataset(fid, "turbulence/u1", datatype(Float32),
    dataspace(n,m,num_of_fields), chunk=(n,m,1), compress=3)
dset_u2 = create_dataset(fid, "turbulence/u2", datatype(Float32),
    dataspace(n,m,num_of_fields), chunk=(n,m,1), compress=3)
dset_temperature = create_dataset(fid, "turbulence/temperature", datatype(Float32),
    dataspace(n,m,num_of_fields), chunk=(n,m,1), compress=3)
dset_q = create_dataset(fid, "turbulence/q", datatype(Float32),
    dataspace(n,m,num_of_fields), chunk=(n,m,1), compress=3)
for i = 1:num_of_fields
    # Read u1
    datatemp = h5read(string("field",i,".h5"),"u1")
    # Save to new file
    dset_u1[:,:,i] = datatemp

    # Read u2
    datatemp = h5read(string("field",i,".h5"),"u2")
    # Save to new file
    dset_u2[:,:,i] = datatemp

    # Read temperature
    datatemp = h5read(string("field",i,".h5"),"temperature")
    # Save to new file
    dset_temperature[:,:,i] = datatemp

    # Read q
    datatemp = h5read(string("field",i,".h5"),"q")
    # Save to new file
    dset_q[:,:,i] = datatemp
end
close(fid)
