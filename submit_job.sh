#!/bin/bash

# Set the number of nodes
#SBATCH --nodes=1

# Set the number of threads per node
#SBATCH --ntasks-per-node=32

# Set max wallclock time
#SBATCH --time=200:00:00

# Set name of job
#SBATCH --job-name=stringEvol

# Mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# Send mail to this address
#SBATCH --mail-user=steven.cotterill@postgrad.manchester.ac.uk

# Set memory requirements
#SBATCH --mem=120GB

# Finish setting up slurm. Now compile and run code

module load gcc7.3.0

#g++ /home/sjcott/Documents/Axions/evolution.cpp /home/sjcott/Documents/Axions/header.cpp -o /home/sjcott/Documents/Axions/Executables/evolution -fopenmp -O3

#/home/sjcott/Documents/Axions/Executables/evolution

g++ /home/sjcott/Documents/Axions/ic_w_evolution.cpp /home/sjcott/Documents/Axions/array.cpp /home/sjcott/Documents/Axions/ic_func.cpp -o /home/sjcott/Documents/Axions/Executables/ic_w_evolution -fopenmp -O3

/home/sjcott/Documents/Axions/Executables/ic_w_evolution


# Now copy output back to telesto

scp /home/sjcott/Documents/Axions/Data/valsPerLoop_loop801.txt sjcott@telesto.jb.man.ac.uk:/mirror2/scratch/sjcott/Documents/C++/Axions/Data/valsPerLoop_loop801.txt


# Finally delete the output file (any slurm... files)

rm slurm-${SLURM_JOB_ID}.out




