#!/bin/bash -l
#SBATCH -A snic2018-8-339
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH -J RNA_assembly_reduced_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user caitlinvigetun@gmail.com
# Load modules
module load bioinfo-tools
module load trinity/2.4.0
module load conda/latest
module load snakemake 
module load samtools/1.8
# Your commands
snakemake
