#!/bin/bash
#
# Minimal working example batch script
#
#SBATCH --partition=short
#SBATCH --job-name=PandemicModelSHP
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=8

python SIRS_grid.py