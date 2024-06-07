#!/bin/bash
#SBATCH -A TSVETANOV-SL3-CPU
#! icelake ccalake
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=00:20:00
#SBATCH -J genfi_cvr_analysis_ca_slurm      # Name of the job
#SBATCH -o genfi_cvr_analysis_ca_slurm.out
#SBATCH -e genfi_cvr_analysis_ca_slurm.err


#! Set the number of subejcts and number of permutations. This will determine the permutation matrix
#! num_subjects=20            # Set this to the number of subejcts in the analysis
num_permutations=10       # Set number of permutations, e.g. 1000
Model="f_rsfa ~ Age + c_Gender" # Set model using Wilkinson Notation (refer to CommonalityAnalysis toolbox for more information)

#! set array number
#SBATCH --array=1-${num_permutations} # Array job to create 100 tasks. number of permutations

#! Modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
#! set up necessary modules
. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
#! module load rhel7/default-ccl               # REQUIRED - loads the basic environment CCLAKE
module load rhel8/default-icl            # REQUIRED - loads the basic environment ICELAKE
module load matlab
#! Path to the MATLAB script. Assuming the currently exectued sh script is in the same directory
script_dir=$(dirname "$(realpath "$0")") # Get the directory of the currently executed script
export PATH=$script_dir:$PATH #! Add MATLAB script path to the PATH environment variable
ca_dir="/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/code/"
export PATH=$ca_dir:$PATH #! Add MATLAB script path to the PATH environment variable

# Set the root directory, which contains subject_info and mask. The results are output in a subfolder in rootDir
rootDir="/rds/project/rds-tbdABxTZvic/genfi_cvr/analysis/vba/ca"
f_data="${rootDir}/subject_info.csv"
#!f_mask="${rootDir}/mask.nii"
f_mask="${rootDir}/SPM12_GM_mask70tpm_91x109x91.nii"

#! Determine number of subjects by counting the number of rows excluding the header in f_data
num_subjects=$(tail -n +2 "$f_data" | wc -l)
echo "Number of subjects: $num_subjects"
echo "Number of permutations: $num_permutations"

#! Initialize output file
date_time=$(date +'%Y%m%d_%H%M') # Get the current date and time in the format YYYYMMDD_HHMM
f_permMatrix="${rootDir}/permutation_matrix_sub${num_subjects}_perm${num_permutations}_${date_time}.txt"

#! Check if the permutation matrix file exists and delete it
if [ -f "$f_permMatrix" ]; then
    echo "File $f_permMatrix exists. Deleting it."
    rm "$f_permMatrix"
fi

#! Generate the permutation matrix with the first column as 1:num_subjects
temp_file=$(mktemp)
> $temp_file

# Add the original permutation as the first column
seq 1 $num_subjects | tr '\n' ' ' | sed 's/ $/\n/' >> $temp_file

#! Generate the remaining permutations
for ((i = 2; i <= num_permutations; i++)); do
    shuf -i 1-$num_subjects -n $num_subjects | tr '\n' ' ' | sed 's/ $/\n/' >> $temp_file
done

#! Move the temporary file to the final output file
mv $temp_file $f_permMatrix

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

#! Extract the column based on SLURM_ARRAY_TASK_ID
#!column_index=${SLURM_ARRAY_TASK_ID:-1}
column_index=$SLURM_ARRAY_TASK_ID
echo "Using column index: $column_index"

matlab -nodisplay -r "ca_vba_run_slurm('$f_data','$Model','$rootDir','$f_permMatrix',$column_index,'$f_mask'); quit"
