#!/bin/bash
#SBATCH --account=myaccount
#SBATCH --job-name=glycine_F
#SBATCH --output=glycine_F.out
#SBATCH --error=glycine_F.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load gaussian16
g16 glycine_F.com
