#!/bin/bash
# Bash script to run an array of jobs on linux machine

# Argument 1 is the number of processors
# Argument 2 is the number of jobs

# Ex: $ ./job_array.sh 16 2

touch number.inp
for i in {1..$2}
do
   echo $i >> number.inp
   mpiexec -n $1 ./turbogen
   rm number.inp
done
