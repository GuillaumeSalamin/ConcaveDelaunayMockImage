#!/bin/bash
#SBATCH --job-name CreateTriangulation
#SBATCH --output /home/astro/ggsalami/ggsalami/TP4b/ConcaveDelaunayMockImage/output/output/out.o%j
#SBATCH --error /home/astro/ggsalami/ggsalami/TP4b/ConcaveDelaunayMockImage/output/error/error.e%j
#SBATCH --mail-user=guillaume.salamin@epfl.ch
#SBATCH --mail-type=ALL

source /home/astro/ggsalami/.bashrc
conda activate pyenv
python CreateTri.py
conda deactivate