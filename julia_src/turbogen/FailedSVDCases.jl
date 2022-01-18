# Check for failed cases and create slurm job submission script to re-run
# the cases that failed

using DelimitedFiles
using HDF5

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
    if ( isfile(string("svd_k",count,".h5")) )
        # Do nothing
    else
        file_exists = 0
        count = count - 1
    end

    # Is file complete?
    if ( file_exists == 1 )
        try
            test_dset2 = h5read(string("svd_k",count,".h5"),"vonkarmanlength/bouyancy")
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

println("Number of datasets checked: ",count)
println("Number of datasets correct: ",num_dsets)

failure = 1
try
    global failed_sims
    failed_sims = readdlm("Failed_Simulations.txt")
    failed_sims = convert(Array{Int64}, failed_sims)
catch
    global failure
    failure = 0
end

if ( failure != 0 )
    failed_sims_string = string("#SBATCH --array=")
    for i ∈ 1:size(failed_sims,1)
        global failed_sims_string
        if ( i != size(failed_sims,1) )
            failed_sims_string = string(failed_sims_string,failed_sims[i],",")
        else
            failed_sims_string = string(failed_sims_string,failed_sims[i],"%","1")
        end
    end
    # Create slurm script
    open("svd_redo.sh", "w") do io
       write(io, "#!/bin/bash")
       write(io, "\n")
       write(io, "#SBATCH --job-name=svd_redo    # Job name")
       write(io, "\n")
       write(io, "#SBATCH --mail-type=END,FAIL    # Only mail on failure or completion")
       write(io, "\n")
       write(io, "#SBATCH --mail-user=acarr1@ufl.edu    # Email address")
       write(io, "\n")
       write(io, "#SBATCH --account=stevenmiller    # Group account")
       write(io, "\n")
       write(io, "#SBATCH --qos=stevenmiller-b    # QOS")
       write(io, "\n")
       write(io, "#SBATCH --ntasks=1    # Number of MPI tasks")
       write(io, "\n")
       write(io, "#SBATCH --cpus-per-task=1    # Number of cpus per task")
       write(io, "\n")
       write(io, "#SBATCH --mem-per-cpu=400mb    # Memory per cpu")
       write(io, "\n")
       write(io, "#SBATCH --time=00-02:00:00    # Time limit")
       write(io, "\n")
       write(io, "#SBATCH --output=svd_redo.log    # Output file")
       write(io, "\n")
       write(io, failed_sims_string)
       write(io, "\n")
       write(io, "\n")
       write(io, "PER_TASK=1")
       write(io, "\n")
       write(io, "START_NUM=\$(( (\$SLURM_ARRAY_TASK_ID - 1) * \$PER_TASK + 1 ))")
       write(io, "\n")
       write(io, "END_NUM=\$(( \$SLURM_ARRAY_TASK_ID * \$PER_TASK ))")
       write(io, "\n")
       write(io, "echo This is task \$SLURM_ARRAY_TASK_ID, which will do runs \$START_NUM to \$END_NUM")
       write(io, "\n")
       write(io, "for (( run=\$START_NUM; run<=END_NUM; run++ )); do")
       write(io, "\n")
       write(io, "echo This is SLURM task \$SLURM_ARRAY_TASK_ID, run number \$run")
       write(io, "\n")
       write(io, "julia \$(pwd)/svd.jl \$run")
       write(io, "\n")
       write(io, "done")
    end
else
    # Do nothing
end
