#!/usr/bin/env bash
# Bash script to run an array of jobs on linux machine

# Argument 1 is the number of processors
# Argument 2 is the number of jobs

# Ex: $ ./job_array.sh 16 2

for i in $(seq 1 $2)
do
  mpiexec -n $1 ./turbogen $i
done
