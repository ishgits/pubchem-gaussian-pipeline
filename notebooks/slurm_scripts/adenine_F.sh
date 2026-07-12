#!/bin/bash
#SBATCH --account=myaccount
#SBATCH --job-name=adenine_F
#SBATCH --output=adenine_F.out
#SBATCH --error=adenine_F.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load gaussian16
g16 adenine_F.com
