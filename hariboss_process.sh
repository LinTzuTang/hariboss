#!/bin/bash
#SBATCH --job-name=hariboss_process    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --partition=bigmem
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=hariboss_20240621_error_revised_process.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python hariboss_process_whole.py

date