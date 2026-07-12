#!/bin/bash
#SBATCH --account=myaccount
#SBATCH --job-name=water_F
#SBATCH --output=water_F.out
#SBATCH --error=water_F.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load gaussian16
g16 water_F.com
