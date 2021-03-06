

using DelimitedFiles

# io = open("waveform_psf.inp", "w")
p = readdlm("waveform_psf.inp")
# close(io)

# Convert to pascals
p_new = p[:,2].*47.88

io = open("F18_waveform.inp", "w")
writedlm(io, [p[:,1] p_new])
close(io)
