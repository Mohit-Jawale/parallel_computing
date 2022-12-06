#!/bin/sh
#
#SBATCH --job-name=dijkstra
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --output=ass_3_20000-%j.out
#SBATCH --time=1000:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjawale@scu.edu
#
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=spread


srun ./dijkstra /WAVE/projects/COEN-319-Fa22/data/pr3/20000.graph 100 output_20000