#!/bin/bash
#SBATCH -A TSVETANOV-SL3-CPU
#! icelake ccalake
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH -J vba_analysis_slurm      # Name of the job
#SBATCH -o vba_analysis_slurm.out
#SBATCH -e vba_analysis_slurm.err
#! set array number
#SBATCH --array=1-200 # Set to match num_permutations below

#! Estimated maximum memory needed (job is force-stopped if exceeded):
#! RAM is allocated in 3420mb blocks, you are charged per block used,
#! and unused fractions of blocks will not be usable by others.
#SBATCH --mem=16G

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
##SBATCH --mail-type=ALL

#! Don't put any #SBATCH directives below this line

# =============================================================================
# CONFIGURATION
# =============================================================================

num_permutations=200       # Must match --array upper bound above

random_seed=12345 # Define a consistent random seed. This is critical to ensure reproducible permutation matrix across nodes/tasks.

column_index=$SLURM_ARRAY_TASK_ID

# Wilkinson notation model: response ~ predictor1 + predictor2 + ...
# Variable names must match column headers in subject_info.xlsx exactly
Model="f_rsfa ~ age + c_sex" # Set model using Wilkinson Notation and variables names in database(refer to CommonalityAnalysis toolbox for more information)


# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
# Modify the environment seen by the application
# (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
# set up necessary modules
. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
#! module load rhel7/default-ccl               # REQUIRED - loads the basic environment CCLAKE
module load rhel8/default-icl            # REQUIRED - loads the basic environment ICELAKE
module load matlab spm

# Paths
script_dir=$(dirname "$(realpath "$0")") # Get the directory of the currently executed script. # Path to the MATLAB script. Assuming the currently exectued sh script is in the same directory
export PATH=$script_dir:$PATH #! Add MATLAB script path to the PATH environment variable

vba_dir="/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/code/"
export PATH=$vba_dir:$PATH #! Add MATLAB script path to the PATH environment variable

spm_dir="/home/kat35/rds/hpc-work/projects/external/mat/spm12_7771/"

# Set the root directory. The results are output in a subfolder in rootDir
outDir="/rds/project/rds-tbdABxTZvic/users/kat35/vba_test/"

f_data="${outDir}data/subject_info.csv"
f_mask="${outDir}data/mask.nii"
#f_mask="${outDir}/data/SPM12_GM_mask70tpm_91x109x91.nii"

echo "Number of permutations: $num_permutations"
echo "Model:                  $Model"

# =============================================================================
# RUN MATLAB VOXEL-BASED ANALYSIS IN MATLAB
# =============================================================================

matlab -nodisplay -nosplash -r \
    "addpath(genpath('${vba_dir}'));addpath(genpath('${spm_dir}')); vba_run_slurm('$f_data','$Model','$outDir',$num_permutations, $random_seed,$column_index,'$f_mask'); quit"
