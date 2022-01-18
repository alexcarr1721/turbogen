#!/bin/bash
#SBATCH --job-name=gather             # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acarr1@ufl.edu    # Where to send mail.  Set this to your email address
#SBATCH --account=stevenmiller        # Account
#SBATCH --qos=stevenmiller            # Quality of service (add -b for burst)
#SBATCH --ntasks=1                    # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1             # Number of cores per MPI task
#SBATCH --mem-per-cpu=8gb             # Memory (i.e. RAM) per processor (for julia, need at least 400mb per cpu)
#SBATCH --time=01-00:00:00            # Wall time limit (days-hrs:min:sec)
#SBATCH --output=gather.log           # Output file

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

julia $(pwd)/gather_svd.jl

echo "Date 		= $(date)"
