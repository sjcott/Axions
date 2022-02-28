#!/bin/bash

# Set the number of tasks
#SBATCH --ntasks=160

# Set the memory limit
#SBATCH --mem=23900MB

# Set max wallclock time
#SBATCH --time=200:00:00

# Set name of job
#SBATCH --job-name=stringEvol

# Mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# Send mail to this address
#SBATCH --mail-user=steven.cotterill@postgrad.manchester.ac.uk

# Finish setting up slurm. Now compile and run code

module load gcc7.3.0

############ Normal evolution compile and run ###################

#g++ /home/sjcott/Documents/Axions/evolution.cpp /home/sjcott/Documents/Axions/header.cpp -o /home/sjcott/Documents/Axions/Executables/evolution -fopenmp -O3

#/home/sjcott/Documents/Axions/Executables/evolution

############ Evolution with initial conditions compile and run ###################

#g++ /home/sjcott/Documents/Axions/ic_w_evolution.cpp /home/sjcott/Documents/Axions/array.cpp /home/sjcott/Documents/Axions/ic_func.cpp -o /home/sjcott/Documents/Axions/Executables/ic_w_evolution -fopenmp -O3

#/home/sjcott/Documents/Axions/Executables/ic_w_evolution

################# MPI evolution compile and run ##################################

mpic++ /home/sjcott/Documents/Axions/mpi_evolution.cpp -o /home/sjcott/Documents/Axions/Executables/mpi_evolution -O3

mpirun -use-hwthread-cpus Executables/mpi_evolution


# Now copy output back to telesto

#scp /home/sjcott/Documents/Axions/Data/valsPerLoop_loop801.txt sjcott@telesto.jb.man.ac.uk:/mirror2/scratch/sjcott/Documents/C++/Axions/Data/valsPerLoop_loop801.txt

scp /home/sjcott/Documents/Axions/GifData/gifStringPosData_G1_PRS_dx0p5_nx1001* sjcott@telesto.jb.man.ac.uk:/mirror2/scratch/sjcott/Documents/C++/Axions/Centaurus/GifData/


# Finally delete the output file (any slurm... files)

rm slurm-${SLURM_JOB_ID}.out




