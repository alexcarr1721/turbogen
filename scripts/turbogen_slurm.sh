#!/bin/bash
#SBATCH --job-name=mpi_run_turbogen       # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acarr1@ufl.edu    # Where to send mail.  Set this to your email address
#SBATCH --account=stevenmiller        # Account
#SBATCH --qos=stevenmiller            # Quality of service (add -b for burst)
#SBATCH --ntasks=16                   # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1             # Number of cores per MPI task 
###SBATCH --nodes=4                     # Maximum number of nodes to be allocated
###SBATCH --ntasks-per-node=32          # Maximum number of tasks on each node
###SBATCH --ntasks-per-socket=16        # Maximum number of tasks on each socket
#SBATCH --mem-per-cpu=800mb           # Memory (i.e. RAM) per processor
#SBATCH --time=00:15:00               # Wall time limit (days-hrs:min:sec)
#SBATCH --output=mpi_run_turbogen.log
#SBATCH --array=1-500%15               # Array range

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module load intel/2020.0.166 openmpi/4.0.3 hdf5/1.12.0

#Set the number of runs that each SLURM task should do
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  srun --mpi=pmix_v3 /blue/stevenmiller/acarr1/software/turbogen/bin/turbogen $run
done
