#!/bin/bash
#SBATCH -A MYACCOUNT-CHANGEME
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:01:00
#SBATCH -J matlab-example
#SBATCH -o matlab_example.out
#SBATCH -e matlab_example.err

module purge
module load rhel7/default-peta4 matlab

nsub=20
nperm=1000

matlab -nodisplay -r "example('output_file'); quit"
