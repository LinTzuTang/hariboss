#!/bin/bash
#SBATCH --job-name=filter_move_results_core    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=96gb                   # Job memory request
#SBATCH --partition=bigmem
#SBATCH --qos=yanjun.li-b
#SBATCH --time=96:00:00             # Time limit hrs:min:sec
#SBATCH --output=filter_move_results_core.log   # Standard output and error log
pwd; hostname; date

module load conda
conda activate hariboss

python filter_move_results.py

date