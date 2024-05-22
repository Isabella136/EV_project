#!/bin/bash

#SBATCH --job-name=ev_analysis
#SBATCH --output=log/ev_analysis.out.%j
#SBATCH --error=log/ev_analysis.err.%j
#SBATCH --time=1-00:00:00
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --directory workdir --unlock
snakemake --profile ./profile/default