#!/bin/bash
# Bash script to run an array of jobs on linux machine

# Argument 1 is the number of processors
# Argument 2 is the number of jobs

# Ex: $ ./job_array.sh 16 2

touch number.inp
counter=0
for i in {1..$2}
do
   let counter++
   mpiexec -n $1 ./turbogen $counter
done
