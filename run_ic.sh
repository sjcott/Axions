#!/bin/bash

# Set the number of tasks
#SBATCH --ntasks=1

# Set name of job
#SBATCH --job-name=ic

module load gcc7.3.0

g++ /home/sjcott/Documents/Axions/initialConditions.cpp /home/sjcott/Documents/Axions/array.cpp -o /home/sjcott/Documents/Axions/Executables/initialConditions -fopenmp -O3

/home/sjcott/Documents/Axions/Executables/initialConditions


