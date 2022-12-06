#!/bin/sh

module switch GCCcore GCC

module load OpenMPI

git pull https://github.com/Mohit-Jawale/parallel_computing.git

make

