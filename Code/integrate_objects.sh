#!/bin/bash

#SBATCH --mail-user=aminerva@princeton.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --time=30:00:00
#SBATCH --mem=240000
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --output='/jukebox/pena/Addie/scRNA-seq/social_vta/logs/%x_%j.log'

module load seurat

Rscript --vanilla integrate_objects.R 