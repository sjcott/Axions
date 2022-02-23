#!/bin/bash
source ~/.bash_functions

mpicmp C++/Axions/ mpi_evolution
if [ $? -ne 0 ]
then
	echo "Compile failed"
else
	mpirun --use-hwthread-cpus Executables/mpi_evolution
	#mpirun Executables/mpi_evolution
fi
