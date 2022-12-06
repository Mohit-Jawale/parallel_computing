#!/bin/sh

module switch GCCcore GCC

module load OpenMPI

git pull https://github.com/Mohit-Jawale/parallel_computing.git

make dijkstra


#
#SBATCH --job-name=matmul
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --output=ass_3-%j.out
#SBATCH --time=1000:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjawale@scu.edu
#
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread



srun ./matmult_omp 1000 1000 1000 20











make clean



