# Post-Processing Code for turbogen

push!(LOAD_PATH, pwd())

using DelimitedFiles
using FFTW
using Interpolations
using HDF5
using Plots
gr()

# Need to load data from each file starting with the first one

# Clear text files
io = open("Successful_Simulations.txt", "w")
close(io)
io = open("Failed_Simulations.txt", "w")
close(io)

# Count datasets
count = 0
file_exists = 1
num_dsets = 0
while file_exists == 1
    global count, file_exists, e1, e2, num_dsets, file_exists
    count = count + 1

    # Does the file exist?
    try
        if ( isfile(string("field",count,".h5")) )
            # Do nothing
        else
            file_exists = 0
        end
    catch
        file_exists = 0
    end

    # Is file complete?
    if ( file_exists == 1 )
        try
            test_dset2 = h5read(string("field",count,".h5"),"microphones/location")
            num_dsets = num_dsets + 1
            open("Successful_Simulations.txt", "a") do io
                writedlm(io, count)
            end
        catch e2
            open("Failed_Simulations.txt", "a") do io
                writedlm(io, count)
            end
        end
    else
        # Do Nothing
    end

end

successful_sims = readdlm("Successful_Simulations.txt")
successful_sims = convert(Array{Int64}, successful_sims)

# Initialize datasets ##########################################################
x₁ = h5read(string("field",successful_sims[1],".h5"), "grid/x1")
x₂ = h5read(string("field",successful_sims[1],".h5"), "grid/x2")
x₃ = h5read(string("field",successful_sims[1],".h5"), "grid/x3")
τ  = h5read(string("field",successful_sims[1],".h5"), "grid/t")
x1half = Int(size(x₁,1)/2)
x2half = Int(size(x₂,1)/2)
x3half = Int(size(x₃,1)/2)
# Contours
ux1_x2x3    = zeros(typeof(x₁[1]), size(x₂,1), size(x₃,1))
ux2_x2x3    = zeros(typeof(x₁[1]), size(x₂,1), size(x₃,1))
ux3_x2x3    = zeros(typeof(x₁[1]), size(x₂,1), size(x₃,1))
T_x2x3      = zeros(typeof(x₁[1]), size(x₂,1), size(x₃,1))
ux1_x2x1    = zeros(typeof(x₁[1]), size(x₂,1), size(x₁,1))
ux2_x2x1    = zeros(typeof(x₁[1]), size(x₂,1), size(x₁,1))
ux3_x2x1    = zeros(typeof(x₁[1]), size(x₂,1), size(x₁,1))
T_x2x1      = zeros(typeof(x₁[1]), size(x₂,1), size(x₁,1))
# Lines
ux1_x1_center   = zeros(typeof(x₁[1]), size(x₁,1))
ux2_x1_center   = zeros(typeof(x₁[1]), size(x₁,1))
ux3_x1_center   = zeros(typeof(x₁[1]), size(x₁,1))
T_x1_center     = zeros(typeof(x₁[1]), size(x₁,1))
ux1_x2_center   = zeros(typeof(x₁[1]), size(x₂,1))
ux2_x2_center   = zeros(typeof(x₁[1]), size(x₂,1))
ux3_x2_center   = zeros(typeof(x₁[1]), size(x₂,1))
T_x2_center     = zeros(typeof(x₁[1]), size(x₂,1))
ux1_x3_center   = zeros(typeof(x₁[1]), size(x₃,1))
ux2_x3_center   = zeros(typeof(x₁[1]), size(x₃,1))
ux3_x3_center   = zeros(typeof(x₁[1]), size(x₃,1))
T_x3_center     = zeros(typeof(x₁[1]), size(x₃,1))

case = "Case1"
fid = h5open(string("TurbogenData",case,".h5"),"w")
create_group(fid, "grid")
g0 = fid["grid"]
g0["tau"] = τ
g0["x1"] = x₁
g0["x2"] = x₂
g0["x3"] = x₃
create_group(fid, "contours_x1_0")
g1 = fid["contours_x1_0"]
dset_ux1_x2x3 = create_dataset(g1, "ux1", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₃,1),num_dsets), chunk=(size(x₂,1),size(x₃,1),1),
    compress=6)
dset_ux2_x2x3 = create_dataset(g1, "ux2", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₃,1),num_dsets), chunk=(size(x₂,1),size(x₃,1),1),
    compress=6)
dset_ux3_x2x3 = create_dataset(g1, "ux3", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₃,1),num_dsets), chunk=(size(x₂,1),size(x₃,1),1),
    compress=6)
dset_T_x2x3 = create_dataset(g1, "T", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₃,1),num_dsets), chunk=(size(x₂,1),size(x₃,1),1),
    compress=6)
create_group(fid, "centerline_x1")
g2 = fid["centerline_x1"]
dset_ux1_x1_center = create_dataset(g2, "ux1", datatype(typeof(x₁[1])),
    dataspace(size(x₁,1),num_dsets), chunk=(size(x₁,1),1),
    compress=6)
dset_ux2_x1_center = create_dataset(g2, "ux2", datatype(typeof(x₁[1])),
    dataspace(size(x₁,1),num_dsets), chunk=(size(x₁,1),1),
    compress=6)
dset_ux3_x1_center = create_dataset(g2, "ux3", datatype(typeof(x₁[1])),
    dataspace(size(x₁,1),num_dsets), chunk=(size(x₁,1),1),
    compress=6)
dset_T_x1_center = create_dataset(g2, "T", datatype(typeof(x₁[1])),
    dataspace(size(x₁,1),num_dsets), chunk=(size(x₁,1),1),
    compress=6)
create_group(fid, "centerline_x2")
gx2 = fid["centerline_x2"]
dset_ux1_x2_center = create_dataset(gx2, "ux1", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),num_dsets), chunk=(size(x₂,1),1),
    compress=6)
dset_ux2_x2_center = create_dataset(gx2, "ux2", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),num_dsets), chunk=(size(x₂,1),1),
    compress=6)
dset_ux3_x2_center = create_dataset(gx2, "ux3", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),num_dsets), chunk=(size(x₂,1),1),
    compress=6)
dset_T_x2_center = create_dataset(gx2, "T", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),num_dsets), chunk=(size(x₂,1),1),
    compress=6)
create_group(fid, "centerline_x3")
gx3 = fid["centerline_x3"]
dset_ux1_x3_center = create_dataset(gx3, "ux1", datatype(typeof(x₁[1])),
    dataspace(size(x₃,1),num_dsets), chunk=(size(x₃,1),1),
    compress=6)
dset_ux2_x3_center = create_dataset(gx3, "ux2", datatype(typeof(x₁[1])),
    dataspace(size(x₃,1),num_dsets), chunk=(size(x₃,1),1),
    compress=6)
dset_ux3_x3_center = create_dataset(gx3, "ux3", datatype(typeof(x₁[1])),
    dataspace(size(x₃,1),num_dsets), chunk=(size(x₃,1),1),
    compress=6)
dset_T_x3_center = create_dataset(gx3, "T", datatype(typeof(x₁[1])),
    dataspace(size(x₃,1),num_dsets), chunk=(size(x₃,1),1),
    compress=6)
create_group(fid, "contours_x3_0")
g3 = fid["contours_x3_0"]
dset_ux1_x2x1 = create_dataset(g3, "ux1", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₁,1),num_dsets), chunk=(size(x₂,1),size(x₁,1),1),
    compress=6)
dset_ux2_x2x1 = create_dataset(g3, "ux2", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₁,1),num_dsets), chunk=(size(x₂,1),size(x₁,1),1),
    compress=6)
dset_ux3_x2x1 = create_dataset(g3, "ux3", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₁,1),num_dsets), chunk=(size(x₂,1),size(x₁,1),1),
    compress=6)
dset_T_x2x1 = create_dataset(g3, "T", datatype(typeof(x₁[1])),
    dataspace(size(x₂,1),size(x₁,1),num_dsets), chunk=(size(x₂,1),size(x₁,1),1),
    compress=6)
# Microphones
# mic_test = h5read(string("field",successful_sims[1],".h5"),"microphones/mics")
# micloc_test = h5read(string("field",successful_sims[1],".h5"),"microphones/location")
# mics    = zeros(typeof(x₁[1]), size(mic_test,1), size(mic_test,2))
# mic_loc = zeros(typeof(x₁[1]), size(micloc_test,1), size(micloc_test,2))
# ASEL    = zeros(typeof(x₁[1]), size(mic_test,2))
# BSEL    = zeros(typeof(x₁[1]), size(mic_test,2))
# CSEL    = zeros(typeof(x₁[1]), size(mic_test,2))
# DSEL    = zeros(typeof(x₁[1]), size(mic_test,2))
# ESEL    = zeros(typeof(x₁[1]), size(mic_test,2))
# PL      = zeros(typeof(x₁[1]), size(mic_test,2))
# ISBAP   = zeros(typeof(x₁[1]), size(mic_test,2))
# create_group(fid, "microphones")
# g4 = fid["microphones"]
# dset_microphones = create_dataset(g4, "mics", datatype(typeof(x₁[1])),
#     dataspace(size(mic_test,1),size(mic_test,2),num_dsets),
#     chunk=(size(mic_test,1),size(mic_test,2),1),
#     compress=6)
# dset_micloc = create_dataset(g4, "mic_loc", datatype(typeof(x₁[1])),
#     dataspace(size(micloc_test,1),size(micloc_test,2),num_dsets),
#     chunk=(size(micloc_test,1),size(micloc_test,2),1),
#     compress=6)

j = 0
for i ∈ successful_sims
    global j, x₁, x₂, x₃, x1half, x2half, x3half, ux1_x2x3, ux2_x2x3,
        ux3_x2x3, T_x2x3, ux1_x1_center, ux2_x1_center,
        ux3_x1_center, T_x1_center, ux1_x2x1, ux2_x2x1, ux3_x2x1, T_x2x1,
        ux1_x2_center, ux2_x2_center, ux3_x2_center, ux1_x3_center,
        ux2_x3_center, ux3_x3_center, T_x2_center, T_x3_center
    j = j + 1

    # Spline interpolation of contours #########################################
    # A3D = h5read(string("field",i,".h5"), "pressure/maximum", (:,:,x1half-2:x1half+3))
    # itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    # range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    # range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    # range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    # sitp = scale(itp, range1, range2, range3 )
    # pmax_x2x3 = sitp(x₂, x₃, (x₁[1]+x₁[size(x₁,1)])/2.0)
    # dset_pmax_x2x3[:,:,j] = pmax_x2x3

    A3D = h5read(string("field",i,".h5"), "turb/ux1", (:,:,x1half-2:x1half+3))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux1_x2x3 = sitp(x₂, x₃, (x₁[1]+x₁[size(x₁,1)])/2.0)
    dset_ux1_x2x3[:,:,j] = ux1_x2x3

    A3D = h5read(string("field",i,".h5"), "turb/ux2", (:,:,x1half-2:x1half+3))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux2_x2x3 = sitp(x₂, x₃, (x₁[1]+x₁[size(x₁,1)])/2.0)
    dset_ux2_x2x3[:,:,j] = ux2_x2x3

    A3D = h5read(string("field",i,".h5"), "turb/ux3", (:,:,x1half-2:x1half+3))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux3_x2x3 = sitp(x₂, x₃, (x₁[1]+x₁[size(x₁,1)])/2.0)
    dset_ux3_x2x3[:,:,j] = ux3_x2x3

    A3D = h5read(string("field",i,".h5"), "turb/temperature", (:,:,x1half-2:x1half+3))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    T_x2x3 = sitp(x₂, x₃, (x₁[1]+x₁[size(x₁,1)])/2.0)
    dset_T_x2x3[:,:,j] = T_x2x3

    # A3D = h5read(string("field",i,".h5"), "pressure/maximum", (:,x3half-2:x3half+3,:))
    # itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    # range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    # range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    # range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    # sitp = scale(itp, range1, range2, range3 )
    # pmax_x2x1 = sitp(x₂, (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    # dset_pmax_x2x1[:,:,j] = pmax_x2x1

    A3D = h5read(string("field",i,".h5"), "turb/ux1", (:,x3half-2:x3half+3,:))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux1_x2x1 = sitp(x₂, (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux1_x2x1[:,:,j] = ux1_x2x1

    A3D = h5read(string("field",i,".h5"), "turb/ux2", (:,x3half-2:x3half+3,:))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux2_x2x1 = sitp(x₂, (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux2_x2x1[:,:,j] = ux2_x2x1

    A3D = h5read(string("field",i,".h5"), "turb/ux3", (:,x3half-2:x3half+3,:))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux3_x2x1 = sitp(x₂, (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux3_x2x1[:,:,j] = ux3_x2x1

    A3D = h5read(string("field",i,".h5"), "turb/temperature", (:,x3half-2:x3half+3,:))
    itp = interpolate(A3D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    T_x2x1 = sitp(x₂, (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_T_x2x1[:,:,j] = T_x2x1
    ############################################################################

    # Lines ####################################################################

    A2D = h5read(string("field",i,".h5"), "turb/ux1",
        (x2half-2:x2half+3,x3half-2:x3half+3,:))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux1_x1_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,
        (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux1_x1_center[:,j] = ux1_x1_center

    A2D = h5read(string("field",i,".h5"), "turb/ux2",
        (x2half-2:x2half+3,x3half-2:x3half+3,:))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux2_x1_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,
        (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux2_x1_center[:,j] = ux2_x1_center

    A2D = h5read(string("field",i,".h5"), "turb/ux3",
        (x2half-2:x2half+3,x3half-2:x3half+3,:))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux3_x1_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,
        (x₃[1]+x₃[size(x₃,1)])/2.0, x₁)
    dset_ux3_x1_center[:,j] = ux3_x1_center

    A2D = h5read(string("field",i,".h5"), "turb/temperature",
        (x2half-2:x2half+3,x3half-2:x3half+3,:))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[1], stop=x₁[size(x₁,1)], length=size(x₁,1))
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    T_x1_center = (sitp((x₂[1] + x₂[size(x₂,1)])/2.0,
        (x₃[1]+x₃[size(x₃,1)])/2.0, x₁).^2)
    dset_T_x1_center[:,j] = T_x1_center

    A2D = h5read(string("field",i,".h5"), "turb/ux1",
        (:,x3half-2:x3half+3,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux1_x2_center = sitp(x₂,
        (x₃[1]+x₃[size(x₃,1)])/2.0, (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux1_x2_center[:,j] = ux1_x2_center

    A2D = h5read(string("field",i,".h5"), "turb/ux2",
        (:,x3half-2:x3half+3,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux2_x2_center = sitp(x₂,
        (x₃[1]+x₃[size(x₃,1)])/2.0, (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux2_x2_center[:,j] = ux2_x2_center

    A2D = h5read(string("field",i,".h5"), "turb/ux3",
        (:,x3half-2:x3half+3,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux3_x2_center = sitp(x₂,
        (x₃[1]+x₃[size(x₃,1)])/2.0, (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux3_x2_center[:,j] = ux3_x2_center

    A2D = h5read(string("field",i,".h5"), "turb/temperature",
        (:,x3half-2:x3half+3,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[1], stop=x₂[size(x₂,1)], length=size(x₂,1))
    range2 = range(x₃[x3half-2], stop=x₃[x3half+3], length=6)
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    T_x2_center = sitp(x₂,
        (x₃[1]+x₃[size(x₃,1)])/2.0, (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_T_x2_center[:,j] = T_x2_center

    A2D = h5read(string("field",i,".h5"), "turb/ux1",
        (x2half-2:x2half+3,:,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux1_x3_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,  x₃,
        (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux1_x3_center[:,j] = ux1_x3_center

    A2D = h5read(string("field",i,".h5"), "turb/ux2",
        (x2half-2:x2half+3,:,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux2_x3_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,  x₃,
        (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux2_x3_center[:,j] = ux2_x3_center

    A2D = h5read(string("field",i,".h5"), "turb/ux3",
        (x2half-2:x2half+3,:,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    ux3_x3_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,  x₃,
        (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_ux3_x3_center[:,j] = ux3_x3_center

    A2D = h5read(string("field",i,".h5"), "turb/temperature",
        (x2half-2:x2half+3,:,x1half-2:x1half+3))
    itp = interpolate(A2D, BSpline(Cubic(Line(OnGrid()))))
    range1 = range(x₂[x2half-2], stop=x₂[x2half+3], length=6)
    range2 = range(x₃[1], stop=x₃[size(x₃,1)], length=size(x₃,1))
    range3 = range(x₁[x1half-2], stop=x₁[x1half+3], length=6)
    sitp = Interpolations.scale(itp, range1, range2, range3 )
    T_x3_center = sitp((x₂[1] + x₂[size(x₂,1)])/2.0,  x₃,
        (x₁[1] + x₁[size(x₁,1)])/2.0 )
    dset_T_x3_center[:,j] = T_x3_center
    ############################################################################
end

close(fid)
